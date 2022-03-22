version development

workflow vamb{

  input{
  ### Inputs to workflow 
  Array[File] FQ1
  Array[File] FQ2
  Array[File] contigs
  Array[String] samples
  File fastxfilter
  }

  
  ### Filter each contig file first
  Array[Pair[String, File]] SampleContigs = zip(samples, contigs) 

  scatter (pair in SampleContigs){
    call filter_assembly{
          input:
            sample = pair.left,
            contigs = pair.right,
            filterscript = fastxfilter
    }
  }
  
  call concatenate_assembly_index{
        input:
          filtered_contigs = filter_assembly.filtered_contigs
  }


  ### Run mapping against the combined contig-catalog
  #Array[Pair[File,File]] Reads = zip(FQ1, FQ2)

  Array[Int] sampleindex = range(length(samples))

  scatter (i in sampleindex){
    File sampleFQ1 = FQ1[i]
    File sampleFQ2 = FQ2[i]
    String sample = samples[i]
  
    call minimap{
        input:
           FQ1 = sampleFQ1,
           FQ2 = sampleFQ2, 
           mmi = concatenate_assembly_index.mmi_index,
           dict = concatenate_assembly_index.contigsdict,
           sample = sample
    }

    call jgidepth{
        input:
           bamfile = minimap.bamfile,
           sample = sample
    }
  }

  File jgi1to3 = jgidepth.jgicontigs[1]

  call vamb{
        input:    
          allcontigs = concatenate_assembly_index.allcontigs,
          jgicounts = jgidepth.jgicounts,
          jgi1to3 = jgi1to3
  }

  output {
    File vambrun = vamb.vamboutput
    File contigs = concatenate_assembly_index.allcontigs
    File contigsindex = concatenate_assembly_index.faifile
  }
}


### Filter by size
task filter_assembly{
  input{
    String sample
    File contigs
    File filterscript
  }
  command{
    python ${filterscript} --i ${contigs} --o ${sample}.assembly.renamed.fna --id ${sample} --min 2000
    pigz ${sample}.assembly.renamed.fna

  }
  output{
    File filtered_contigs = "${sample}.assembly.renamed.fna.gz"
  }

  runtime{
      docker: "gcr.io/microbiome-xavier/metagenomicstools:032819"
      disks: "local-disk 40 HDD"
      cpu : "1"
      preemptible: 2
  }

}

### Concatenate
task concatenate_assembly_index{
  input{
    Array[File] filtered_contigs
  }
  command{
      cat ~{sep=' ' filtered_contigs} > allcontigs.fna.gz 

      echo "staring minimap2 index"
      minimap2 -I 12G -d contigs.flt.mmi allcontigs.fna.gz  2> index_log.txt

      echo "Creating contig dictionary with samtools"
      samtools dict allcontigs.fna.gz | cut -f1-3 > contigs.flt.dict 2> log.txt
      ls
      pigz -dc allcontigs.fna.gz > allcontigs.fna
      samtools faidx allcontigs.fna
  }
  output{
    File allcontigs = "allcontigs.fna.gz"
    File contigsdict = "contigs.flt.dict"
    File mmi_index = "contigs.flt.mmi"
    File faifile = "allcontigs.fna.fai"
  }

  runtime{
      docker: "nanozoo/minimap2:latest"
      disks: "local-disk 400 HDD"
      memory: "128GB"
      cpu : "1"
      preemptible: 1
  }

}

task minimap{
  input{
    File FQ1
    File FQ2
    File mmi
    File dict
    String sample
  }
  command{
      minimap2 -t 30 -ax sr ${mmi} ${FQ1} ${FQ2} | grep -v "^@" | cat ${dict} - | samtools view -F 3584 -b - > ${sample}.bam 2> mapping_log.txt
  }
  output{
    File bamfile = "${sample}.bam"

  }
  runtime{
      docker: "nanozoo/minimap2:latest"
      disks: "local-disk 200 HDD"
      memory: "40GB"
      cpu : "30"
      preemptible: 2
  }

}

task jgidepth{
    input{
    File bamfile
    String sample
    }
    command{
        samtools sort ${bamfile} -T tmp.${sample} --threads 4 -m 5G -o ${sample}.sort.bam
        jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth ${sample}.jgidepth ${sample}.sort.bam
        cut -f1-3 ${sample}.jgidepth > ${sample}.jgidepth.1to3
        cut -f1-3 --complement ${sample}.jgidepth > ${sample}.jgidepth.counts
    }
    output{
      File jgidepth = "${sample}.jgidepth"
      File jgicounts = "${sample}.jgidepth.counts"
      File jgicontigs = "${sample}.jgidepth.1to3"
    }
    runtime{
        docker: "biojokke/vamb:latest"
        disks: "local-disk 150 HDD"
        memory: "30GB"
        cpu : "4"
        preemptible: 1
    }
}

task vamb{
  input{
    File allcontigs
    Array[File] jgicounts
    File jgi1to3

  }
  command{

         ### assembly depth matrix 
         paste ${jgi1to3} ${sep=' ' jgicounts} > jgi.abundance.matrix

         vamb -o :: --minfasta 800000 --outdir vamb --fasta ${allcontigs} --jgi jgi.abundance.matrix -n 512 512 -l 32 -b 200

         ### bins can be found in vamb/bins
         gzip jgi.abundance.matrix
         gzip vamb/bins/* 
         mv jgi.abundance.matrix.gz vamb
         tar -zcvf vamb.tar.gz vamb
  }
  output{
    File vamboutput = "vamb.tar.gz"
  }
  runtime{
      docker: "biojokke/vamb:latest"
      disks: "local-disk 500 HDD"
      memory: "120GB"
      cpu : "32"
      preemptible: 1
  }

}

