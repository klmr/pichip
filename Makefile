SHELL := $(shell which bash)
.DELETE_ON_ERROR:

keep = $(foreach i,$2,$(if $(findstring $1,$i),$i))
filter_out = $(foreach i,$2,$(if $(findstring $1,$i),,$i))

faidx = raw/reference/ce10.chrom.sizes

bam = $(shell find ./raw/mapped/SX*-$1/ -name '*.bam' -depth 1 -print)

mkdir = mkdir -p $(dir $@)

wt_files = $(call bam,wt)
hrde1_files = $(call bam,hrde1)
emb4_files = $(call bam,emb4)

raw_files = ${wt_files} ${hrde1_files} ${emb4_files}

cov_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.gc})
bedgraph_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.bg})
rel_bg_files = ${bedgraph_files:.bg=.rel_bg}

data/coverage/%.gc: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -ibam $< -g $(lastword $^) > $@

data/coverage/%.bg: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -bg -ibam $< -g $(lastword $^) > $@

data/coverage/%.rel_bg: data/coverage/%.bg raw/mapped/%.bam
	@# Normalize by per-base coverage scaling factor
	awk -F $$'\t' ' \
		BEGIN { \
			covered_bases = '$$(samtools view -F4 '$(lastword $^)' | awk -F $$'\t' '{covered_bases += length($$10)} END {print covered_bases}')'; \
			genome_length = '$$(awk '{x += $$2} END {print x}' ${faidx})'; \
			sf = covered_bases / genome_length \
		} \
		{ \
			$$4 = $$4 / sf; \
			print $$0 \
		}' $< \
	> $@

.SECONDEXPANSION:

replicates = $(shell find ./data/coverage/$1/ -name '*.rel_bg' -print)
chip_replicates = $(call filter_out,input,$(call replicates,$1))
input_replicates = $(call keep,input,$(call replicates,$1))

data/coverage/%/chip_merged.bg: $$(call chip_replicates,$$*)
	./scripts/mean_bedgraph $^ > $@

data/coverage/%/input_merged.bg: $$(call input_replicates,$$*)
	./scripts/mean_bedgraph $^ > $@

data/coverage/%.bw: data/coverage/%.bg ${faidx}
	tmpfile="$$(mktemp)"; \
	LC_COLLATE=C sort -k1,1 -k2,2n $< > "$$tmpfile"; \
	trap "rm -f $$tmpfile" EXIT; \
	bedGraphToBigWig "$$tmpfile" $(lastword $^) $@
