SHELL := $(shell which bash)

faidx = raw/reference/ce10.chrom.sizes

bam = $(shell find ./raw/mapped/SX*-$1/ -name '*.bam' -depth 1 -print)

mkdir = mkdir -p '$(dir $@)'

wt_files = $(call bam,wt)
hrde1_files = $(call bam,hrde1)
emb4_files = $(call bam,emb4)

raw_files = ${wt_files} ${hrde1_files} ${emb4_files}

cov_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.gc})
bedgraph_files = $(subst raw/mapped/,data/coverage/,${raw_files:.bam=.bg})
rel_bg_files = ${bedgraph_files:.bg=.rel_bg}

data/coverage/%.gc: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -ibam '$<' -g '$(lastword $^)' > '$@'

data/coverage/%.bg: raw/mapped/%.bam raw/mapped/%.bam.bai ${faidx}
	@$(mkdir)
	bedtools genomecov -bg -ibam '$<' -g '$(lastword $^)' > '$@'

data/coverage/%.rel_bg: data/coverage/%.bg raw/mapped/%.bam
	awk ' \
		BEGIN { \
			sf = '$$(samtools view -F4 -c '$(lastword $^)')' / \
				'$$(awk '{ x += $$2 } END { print x }' ${faidx})' \
		} \
		{ \
			$$4 = $$4 / sf; \
			print $$0 \
		}' '$<' \
	> '$@'
