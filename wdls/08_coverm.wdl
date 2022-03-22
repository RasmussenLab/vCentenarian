workflow coverm_map{
	File FQ1
    File FQ2
	File genomes
	String sample


	call map_coverm {
		input:
			fileR1=FQ1,
			fileR2=FQ2,
			sample=sample,
			genomes=genomes
	}

    call collect{
        input:
            covermfiles = map_coverm.fileCoverm
    }

	output{
		File covermfile = collect.coverm_files
	}
}

task map_coverm {
	File fileR1
	File fileR2
	String sample
	File genomes

	command {
        coverm contig --contig-end-exclusion 0 -m count rpkm covered_fraction length --min-read-percent-identity 95 -1 ${fileR1} -2 ${fileR2} \
        --reference ${genomes} \
        --threads 8 \
        --output-file ${sample}.coverm.txt
        
        cut -f1 ${sample}.coverm.txt > gene.names.txt
        cut -f1,2 ${sample}.coverm.txt > ${sample}.coverm.count.txt
	}

	output {
		File fileCount = "${sample}.coverm.count.txt"
        File fileCoverm = "${sample}.coverm.txt"
        File fileGenenames = "gene.names.txt"
	}
	runtime {
		docker: "quay.io/biocontainers/coverm:0.6.1--ha1d52fc_1"
		cpu: 8
  		memory: "32GB"
  		preemptible: 2
  		bootDiskSizeGb: 50
  		disks: "local-disk 300 HDD"
		maxRetries: 2
	}
}

task collect{
  input{
    Array[File] covermfiles
  }
  command <<<

        mkdir coverm
        for f in ~{sep=' ' covermfiles}; do
        mv ${f} coverm
        done
        tar -czvf covermfiles.tar.gz coverm 

  >>>
  output{
    File coverm_files = "covermfiles.tar.gz"
  }
  runtime{
      docker: "ubuntu:16.04"
      disks: "local-disk 150 HDD"
      cpu : "1"
      preemptible: 1
  }

}