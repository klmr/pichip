SHELL := $(shell which bash)
.DELETE_ON_ERROR:

keep = $(foreach i,$2,$(if $(findstring $1,$i),$i))
filter_out = $(foreach i,$2,$(if $(findstring $1,$i),,$i))

faidx = raw/reference/ce10.chrom.sizes

bam = $(shell find ./raw/mapped/$1/ -name '*.bam' -depth 1 -print)

mkdir = mkdir -p $(dir $@)

conditions = $(notdir $(shell find ./raw/mapped -name 'SX*' -depth 1 -print))
raw_files = $(foreach c,${conditions},$(call bam,$c))

cov_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.gc})
raw_bedgraph_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.bedgraph})
rel_bg_files = ${raw_bedgraph_files:.bedgraph=.rel_bg}

data/coverage/%.gc: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -ibam $< -g $(lastword $^) > $@

data/coverage/%.bedgraph: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -bg -ibam $< -g $(lastword $^) > $@

data/coverage/%.rel_bg: data/coverage/%.bedgraph raw/mapped/%.bam
	@# Normalize by per-base coverage scaling factor
	awk -F $$'\t' -v OFS=$$'\t' ' \
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

replicates = $(call keep,$1,${rel_bg_files})
chip_replicates = $(call filter_out,input,$(call replicates,$1))
input_replicates = $(call keep,input,$(call replicates,$1))

chip_stem = $(addsuffix _chip_merged,${conditions})
input_stem = $(addsuffix _input_merged,${conditions})
all_stem = ${chip_stem} ${input_stem}

bedgraph_files = $(addprefix data/coverage/bedgraph/,$(addsuffix .bedgraph,${all_stem}))
.PRECIOUS: ${bedgraph_files}

.PHONY: bedgraph-tracks
## BedGraph tracks for the merged, normalized ChIP and input libraries
bedgraph-tracks: ${bedgraph_files}

data/coverage/bedgraph/%_chip_merged.bedgraph: $$(call chip_replicates,$$*)
	@$(mkdir)
	./scripts/mean_bedgraph $^ > $@

data/coverage/bedgraph/%_input_merged.bedgraph: $$(call input_replicates,$$*)
	@$(mkdir)
	./scripts/mean_bedgraph $^ > $@

bigwig_files = $(addprefix data/coverage/bigwig/,$(addsuffix .bw,${all_stem}))
.PRECIOUS: ${bigwig_files}

.PHONY: bigwig-tracks
## BigWig tracks for the merged, normalized ChIP and input libraries
bigwig-tracks: ${bigwig_files}

data/coverage/bigwig/%.bw: data/coverage/bedgraph/%.bedgraph ${faidx}
	@$(mkdir)
	tmpfile="$$(mktemp)"; \
	LC_COLLATE=C sort -k1,1 -k2,2n $< > "$$tmpfile"; \
	trap "rm -f $$tmpfile" EXIT; \
	bedGraphToBigWig "$$tmpfile" $(lastword $^) $@

bam_files = $(call keep,$1,${raw_files})
chip_bam_files = $(call filter_out,input,$(call bam_files,$1))
input_bam_files = $(call keep,input,$(call bam_files,$1))

data/peaks/%_peaks.bw: $$(call bam_files,$$(dir $$*)) ${faidx}
	macs2 callpeak --broad --treatment $(call chip_bam_files,$(dir $*)) --control $(call input_bam_files,$(dir $*)) --name $(notdir $*) --format BAM --outdir $(dir $@) --bdg
	tmpfile="$$(mktemp)"; \
	LC_COLLATE=C sort -k1,1 -k2,2n $(@:_peaks.bw=_treat_pileup.bdg) > "$$tmpfile"; \
	trap "rm -f $$tmpfile" EXIT; \
	bedGraphToBigWig "$$tmpfile" $(lastword $^) $@

macs2_peaks = $(foreach c,${conditions},data/peaks/$c/$c_peaks.bw)

.PHONY: call-peaks
## Call ChIP peaks using macs2
call-peaks: ${macs2_peaks}

# FIXME: Narrow peaks used for subsequent analysis, but above command created broad peaks.
macs2_narrow_peaks = $(foreach c,${conditions},data/peaks-narrow/$c/$c_peaks.bw)

peak_prefixes = $(subst _peaks.bw,,$(call filter_out,wt,${macs2_narrow_peaks}))
peak_treatment = $(call keep,$1,$(addsuffix _treat_pileup.bdg,${peak_prefixes}))
peak_control = $(call keep,$1,$(addsuffix _control_lambda.bdg,${peak_prefixes}))
all_peak_files = $(call peak_treatment,$1) $(call peak_control,$1)

get_eff_seq_depth = $(shell grep 'tags after filtering in \(treatment\|control\)' ${1:_treat_pileup.bdg=_peaks.xls} | grep -o '[[:digit:]]\+' | sort -n | head -1)

data/differential-peaks/%_common.bed: $$(call all_peak_files,$$(dir $$*))
	@$(mkdir)
	$(eval d1 := $(call get_eff_seq_depth,./data/peaks-narrow/SX1316-wt/SX1361-wt_treat_pileup.bdg))
	$(eval d2 := $(call get_eff_seq_depth,$<))
	macs2 bdgdiff \
		--t1 ./data/peaks-narrow/SX1316-wt/SX1361-wt_treat_pileup.bdg \
		--c1 ./data/peaks-narrow/SX1316-wt/SX1361-wt_control_lambda.bdg \
		--t2 $(firstword $+) --c2 $(lastword $+) \
		--d1 ${d1} --d2 ${d2} -g 100 -l 500 --outdir $(dir $@) \
		-o $(notdir ${@:common.bed=down.bed} ${@:common.bed=up.bed} $@)

differential_peaks = $(foreach c,$(call filter_out,wt,${conditions}),data/differential-peaks/$c/diff_wt-vs-$c_common.bed)

.PHONY: differential-peaks
## Call differential ChIP-seq peaks
differential-peaks: ${differential_peaks}

normalized_bigwig_files = $(addprefix data/coverage/bigwig/,$(addsuffix _merged.bw,${conditions}))

.PHONY: normalized-bigwig-tracks
## BigWig tracks for the merged ChIP libraries, normalized to input
normalized-bigwig-tracks: ${normalized_bigwig_files}

normalized_bedgraph_files = $(addprefix data/coverage/bedgraph/,$(addsuffix _merged.bedgraph,${conditions}))
.PRECIOUS: ${normalized_bedgraph_files}

data/coverage/bedgraph/%_merged.bedgraph: data/coverage/bedgraph/%_chip_merged.bedgraph data/coverage/bedgraph/%_input_merged.bedgraph
	./scripts/normalize_to_input --faidx ${faidx} $+ > $@

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
