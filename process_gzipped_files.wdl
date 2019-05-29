workflow EnrichmentAnalysisWorkflow {
    File geneSampleCoverageControlGz
    File geneSampleCoverageCaseGz
    File intervalSampleCoverageControlGz
    File intervalSampleCoverageCaseGz

	String sampleSetId
    Int intervalsPerBatch=2000
    Int genesPerBatch=500
    String dockerVersion="11.0"

    call PrepareEnrichment_Task as PrepareEnrichment_GeneControl {
        input: sampleCoverageGz=geneSampleCoverageControlGz,
        	   dockerVersion=dockerVersion,
        	   intervalsPerBatch=genesPerBatch
    }

    call PrepareEnrichment_Task as PrepareEnrichment_GeneCase {
        input: sampleCoverageGz=geneSampleCoverageCaseGz,
        	   dockerVersion=dockerVersion,
        	   intervalsPerBatch=genesPerBatch
    }

    call PrepareEnrichment_Task as PrepareEnrichment_IntervalControl {
        input: sampleCoverageGz=intervalSampleCoverageControlGz,
        	   dockerVersion=dockerVersion,
        	   intervalsPerBatch=intervalsPerBatch
    }

    call PrepareEnrichment_Task as PrepareEnrichment_IntervalCase {
        input: sampleCoverageGz=intervalSampleCoverageCaseGz,
        	   dockerVersion=dockerVersion,
        	   intervalsPerBatch=intervalsPerBatch
    }

    scatter (i in range(length(PrepareEnrichment_GeneControl.shards))) {
        call Comparison_Task as Comparison_GeneTask {
            input: controlShard=PrepareEnrichment_GeneControl.shards[i],
            	   caseShard=PrepareEnrichment_GeneCase.shards[i],
                   geneOrInterval="gene",
            	   dockerVersion=dockerVersion
        }
    }

    scatter (i in range(length(PrepareEnrichment_IntervalControl.shards))) {
        call Comparison_Task as Comparison_IntervalTask {
            input: controlShard=PrepareEnrichment_IntervalControl.shards[i],
            	   caseShard=PrepareEnrichment_IntervalCase.shards[i],
                   geneOrInterval="interval",
            	   dockerVersion=dockerVersion
        }
    }

    call GatherShards_Task as GatherShards_GeneTask {
    	input: dockerVersion=dockerVersion,
               shards=Comparison_GeneTask.result,
               sampleSetId=sampleSetId
	}

    call GatherShards_Task as GatherShards_IntervalTask {
    	input: dockerVersion=dockerVersion,
               shards=Comparison_IntervalTask.result,
               sampleSetId=sampleSetId
	}

    output {
    	GatherShards_GeneTask.combinedIntervals
        GatherShards_IntervalTask.combinedIntervals
	}
}

task GatherShards_Task {
    Array[File] shards
    File firstShard = shards[0]
	String sampleSetId
    String dockerVersion

    command <<<
        # Add the column headers into the final file
        echo "cat ${firstShard} | head -n1 > '${sampleSetId}.intervals_combined.txt'"
        cat ${firstShard} | head -n1 > "${sampleSetId}.intervals_combined.txt"

        # Add content from each of the shards (minus the headers) into the final file
        for f in ${sep = ' ' shards}
        do
        	cat $f | tail -n +2 >> "${sampleSetId}.intervals_combined.txt"
        done

        echo "Number of lines in final concatentated file:"
        echo "cat ${sampleSetId}.intervals_combined.txt | wc -l"
        cat ${sampleSetId}.intervals_combined.txt | wc -l
    >>>

    runtime {
        docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
    }

    output {
        File combinedIntervals="${sampleSetId}.intervals_combined.txt"
    }
}


task Comparison_Task {
	String dockerVersion
    File controlShard
    File caseShard
    String geneOrInterval

    command <<<
        mkdir outputs
       	python /differential_coverage_analysis.py --fractions_case ${caseShard} --fractions_control ${controlShard} --output_folder outputs --gene_or_interval ${geneOrInterval}

        echo "ls -lh outputs"
        ls -lh outputs

        echo "mv outputs/* /cromwell_root/"
        mv outputs/* /cromwell_root/

        echo "ls -lh /cromwell_root/"
        ls -lh /cromwell_root/
    >>>

    output {
    	File result="${geneOrInterval}_dc_report.tsv"
    }

    runtime {
        docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
        preemptible: 3
    }
}


task PrepareEnrichment_Task {
    File sampleCoverageGz
    Int intervalsPerBatch
	String dockerVersion
    Int prepare_disk_gb

    command <<<
        set -xeuo pipefail
        function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 30;
       	done &

    	# Unzip the bzipped file
        echo "zcat ${sampleCoverageGz} > coverage.txt"
        zcat ${sampleCoverageGz} > coverage.txt

        echo "ls -lh coverage.txt"
        ls -lh coverage.txt

        # Split the merged file into # 20,000 / intervalsPerBatch files (exclude the header)
        tail -n +2 coverage.txt | split -l ${intervalsPerBatch} - split_

        echo "Word counts after split..."
        for f in $(ls split_*); do echo $f; cat $f | wc -l;  done

        # Add the header back in for each of sharded files
        echo "Adding the header back in...."
        for file in split_*
        do
            head -n 1 coverage.txt > tmp_file
            cat $file >> tmp_file
            mv -f tmp_file $file
        done

        echo "Word counts after adding headers in..."
        for f in $(ls split_*); do echo $f; cat $f | wc -l;  done

        # Keep track of the column headers for use later in this workflow
        head -n 1 coverage.txt > column_headers.txt

        echo "Number of lines in original file:"
        echo "cat coverage.txt | wc -l"
        cat coverage.txt | wc -l
    >>>

    output {
        Array[File] shards=glob('split_*')
        File columnHeaders="column_headers.txt"
    }

    runtime {
        docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
		preemptible: 3
        disks: "local-disk ${prepare_disk_gb} HDD"
    }
}