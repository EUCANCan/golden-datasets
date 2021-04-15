.PHONY: all clean download bam2fastq compress check
.ONESHELL:
.PRECIOUS: %.local.hashes

SHELL   = /bin/bash
WGET    = wget
CONDA   = conda
MKDIR_P = mkdir -p
RM_RF   = rm -rf

TEST               := false

CONDA_BASE          = $(shell conda info --base)
CONDA_ENV           = $(addsuffix /env, $(realpath .))
CONDA_URL           = https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
CONDA_ACTIVATE      = source $(CONDA_BASE)/etc/profile.d/conda.sh ; conda activate ; conda activate

SYNAPSE_SESSION     = $(HOME)/.synapseSession
BEDTOOLS_URL        = https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary

DATA_DIR            = data
DEFAULT_DIR         = $(DATA_DIR) # Change this line if you want to change the default directory

# INSTALL_INPUT	:= $(realpath $(shell read -p "Enter path to installation directory [build]:"$$'\n' && echo $$REPLY))

ifeq ($(TEST),true)
DATASETS_TSV        := $(DATA_DIR)/datasets.test.tsv
INSTALL_DIR         ?= $(or $(INSTALL_INPUT), $(strip $(DEFAULT_DIR)))/test
else
DATASETS_TSV        := $(DATA_DIR)/datasets.tsv
INSTALL_DIR         ?= $(or $(INSTALL_INPUT), $(strip $(DEFAULT_DIR)))
endif

#TODO: Dynamic set doesn't work correctly

EXTRACT_TSV_CMD       := awk 'BEGIN {ext=""} NR<=1 {next} {if($$5 !~ /null/) ext="gsutil"; else if ($$4 !~ /null/) ext="synid"; else ext="url"; printf "%s/%s.%s\n", $$1, $$2, ext}'

DLIDS_FILES           := $(addprefix $(INSTALL_DIR)/, $(shell $(EXTRACT_TSV_CMD) $(DATASETS_TSV)))

URL_FILES             := $(filter %.url, $(DLIDS_FILES))
SYNIDS_FILES          := $(filter %.synid, $(DLIDS_FILES))
GSUTIL_FILES          := $(filter %.gsutil, $(DLIDS_FILES))

DATASETS_FROM_URL_FILES         :=  $(URL_FILES:.url=)
DATASETS_FROM_GSUTIL_FILES      :=  $(GSUTIL_FILES:.gsutil=)
DATASETS_FROM_SYNID_FILES       :=  $(SYNIDS_FILES:.synid=)

DATASETS_FILES        := $(DATASETS_FROM_URL_FILES) $(DATASETS_FROM_SYNID_FILES) $(DATASETS_FROM_GSUTIL_FILES)

DATASETS_DIRS         := $(shell echo "$(dir $(DATASETS_FILES))" | tr ' ' '\n' | uniq)

BAM_FILES             := $(shell find $(INSTALL_DIR)/ -type f -name '*.bam')
BAM2FASTQ_FILES       := $(BAM_FILES:.bam=.r1.fastq) $(BAM_FILES:.bam=.r2.fastq)
FASTQ2GZ_FILES        := $(BAM2FASTQ_FILES:.fastq=.fastq.gz)

CHECKSUM_FILES           := $(GSUTIL_FILES:.gsutil=.checksum) $(SYNIDS_FILES:.synid=.checksum) $(URL_FILES:.url=.checksum)
LOCAL_HASHES_FILES       := $(CHECKSUM_FILES:.checksum=.local.hashes)

ifeq (, $(shell which $(CONDA)))
$(error "No $(CONDA) in $PATH, check miniconda website to install it ($(CONDA_URL))")
endif

###############################################################################
#                                  GOALS
###############################################################################

all: download

download: $(DATASETS_FILES)

bam2fastq: $(BAM2FASTQ_FILES)

compress: $(FASTQ2GZ_FILES)

check: $(CHECKSUM_FILES)

clean:
	$(RM_RF) $(LOCAL_HASHES_FILES)


###############################################################################
#                                 CONDA
###############################################################################

# Create the environment, install and initialize all the tools required
$(CONDA_ENV): environment.yml
	rm -rf $@ $(SYNAPSE_SESSION)
	$(CONDA) env create -p $@ --file $<


###############################################################################
#                                 DOWNLOAD
###############################################################################
.INTERMEDIATE: $(DLIDS_FILES)

%/:
	$(MKDIR_P) $(@D)

%.synid: $(DATASETS_TSV) | $(DATASETS_DIRS)
	@out=$$(awk '$$2 ~ /^$(notdir $*)$$/ {if ($$4 != "null")print $$4; else exit 1}' $<) && echo $$out > $@ || echo "e[33mSynapse ID not found for $*\e[0m"

%.gsutil: $(DATASETS_TSV) | $(DATASETS_DIRS)
	@out=$$(awk '$$2 ~ /^$(notdir $*)$$/ {if ($$5 != "null")print $$5; else exit 1}' $<) && echo $$out > $@  || echo "e[33mGoogle storage URI not found for $*\e[0m"

%.url: $(DATASETS_TSV) | $(DATASETS_DIRS)
	@out=$$(awk '$$2 ~ /^$(notdir $*)$$/ {if ($$3 != "null")print $$3; else exit 1}' $<) && echo $$out > $@  || echo "e[33mURL not found for $*\e[0m"

# If synapse credentials not saved before, login and save them with the OS keyring tool
$(SYNAPSE_SESSION): | $(CONDA_ENV)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	synapse login --rememberMe

# Download file from synapse if there is an URL in the tsv file
# TODO: use rsync instead
$(DATASETS_FROM_URL_FILES): %: | %.url
	$(WGET) -O $@ $(shell cat $(word 1, $|))

# Download file from synapse if there is a synapse ID in the tsv file
$(DATASETS_FROM_SYNID_FILES): %: | %.synid $(SYNAPSE_SESSION)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	syn_file=$$(synapse get --downloadLocation $(@D) $(shell cat $(word 1, $|))  | tee /dev/tty | awk 'BEGIN {FS=" "} /Downloaded file/ {print $$3}') && \
	mv $(@D)/$${syn_file} $@

$(DATASETS_FROM_GSUTIL_FILES): %: | %.gsutil $(CONDA_ENV)
	gsutil cp  $(shell cat $(word 1, $|)) $(@D)


###############################################################################
#                                 bam2fastq
###############################################################################

# bedtools bamtofastq -i $< -fq $(word 1, $@)
%.r1.fastq %.r2.fastq: %.bam | $(CONDA_ENV)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	gatk SamToFastq --INPUT=$< --FASTQ=$*.r1.fastq --SECOND_END_FASTQ=$*.r2.fastq

%.fastq.gz: %.fastq | $(CONDA_ENV)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	pigz $<


###############################################################################
#                                 check
###############################################################################
.INTERMEDIATE: $(CHECKSUM_FILES)

# TODO: use sha256sum
%.ori.hashes: %.synid | $(SYNAPSE_SESSION)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	synapse show $(shell cat $(word 1, $<)) | awk -F '=' '/md5/{print $$2}' > $@

# Get file Hash from TSV or trhough an API like gsutil + URI
%.ori.hashes: %.url | $(CONDA_ENV)
	@echo -e '\e[33mMD5 checksum from url not supported yet (will be added as a new column in another release).\e[0m'
	touch $@

%.ori.hashes: %.gsutil | $(CONDA_ENV)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	gsutil hash $(shell cat $(word 1, $<)) | awk '/crc32c/{print $$3}' > $@

# For Gsutil generate crc32c/md5 hash in hexadecimal format
%.local.hashes: %.gsutil | $(CONDA_ENV)
	@$(CONDA_ACTIVATE) $(CONDA_ENV)
	gsutil hash $* | awk '/crc32c/{print $$3}' > $@

%.local.hashes: %.synid | $(CONDA_ENV)
	md5sum $* | awk '{print $$1}' > $@

%.checksum: %.ori.hashes %.local.hashes | $(CONDA_ENV)
	@diff $^ && echo -e "\e[32mHashes are identical for $@ \e[0m" && touch $@ || echo -e "\e[31mHashes are not the same for $@ !\e[0m"
