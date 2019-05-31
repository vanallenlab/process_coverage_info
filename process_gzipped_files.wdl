workflow process_coverage_info {
	File sampleSummaryFilesGzipped
    File intervalSummaryFilesGzipped
    File geneSummaryFilesGzipped
    Int diskSpace
    String memory
    String label
    Int cutoff=25
    Int preemptible
    Int bootSizeDiskGb=100
    String dockerVersion="28.0"
    Int threads=10

    call ProcessCoverageInfoSampleTask {
    	input:
        	sampleSummaryFilesGzipped=sampleSummaryFilesGzipped,
        	diskSpace=diskSpace,
            memory=memory,
            label=label,
            cutoff=cutoff,
            bootSizeDiskGb=bootSizeDiskGb,
            preemptible=preemptible,
            dockerVersion=dockerVersion
    }

    call ProcessCoverageInfoGeneTask {
    	input:
            geneSummaryFilesGzipped=geneSummaryFilesGzipped,
        	diskSpace=diskSpace,
            memory=memory,
            label=label,
            cutoff=cutoff,
            bootSizeDiskGb=bootSizeDiskGb,
            preemptible=preemptible,
            dockerVersion=dockerVersion,
            threads=threads
    }

    call ProcessCoverageInfoIntervalTask {
    	input:
            intervalSummaryFilesGzipped=intervalSummaryFilesGzipped,
        	diskSpace=diskSpace,
            memory=memory,
            label=label,
            cutoff=cutoff,
            bootSizeDiskGb=bootSizeDiskGb,
            preemptible=preemptible,
            dockerVersion=dockerVersion,
            threads=threads
    }
}

task ProcessCoverageInfoSampleTask {
	File sampleSummaryFilesGzipped
    Int diskSpace
    String memory
    String label
    Int cutoff
    Int preemptible
    String bootSizeDiskGb
    String dockerVersion

    command <<<
    	# log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

        tar -xvzf ${sampleSummaryFilesGzipped}

        echo "ls -lh"
        ls -lh

        echo "mkdir output"
        mkdir output

        echo "ls -lh sample_summary_files | wc"
        ls -lh sample_summary_files | wc

        # Run Python splice junction discovery script
        echo "python /process_coverage_info.py --label ${label} --sample sample_summary_files --cutoff ${cutoff} --output /cromwell_root/"
        python /process_coverage_info.py --label ${label} --sample sample_summary_files --cutoff ${cutoff} --output /cromwell_root/

        echo "Output files in /cromwell_root/:"
        ls -lh /cromwell_root/

    >>>

    output {
		File mean_coverage_by_sample="${label}_mean_coverage_by_sample_id.tsv"
    }

    runtime {
      	docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
        memory: "${memory}"
        disks: "local-disk " + diskSpace + " HDD"
        bootDiskSizeGb: "${bootSizeDiskGb}"
        preemptible: preemptible
    }
}

task ProcessCoverageInfoGeneTask {
    File geneSummaryFilesGzipped
    Int diskSpace
    String memory
    String label
    Int cutoff
    Int preemptible
    String bootSizeDiskGb
    String dockerVersion
	Int threads

    command <<<
    	# log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

        tar -xvzf ${geneSummaryFilesGzipped}

        echo "ls -lh"
        ls -lh

        echo "mkdir output"
        mkdir output

        echo "ls -lh gene_summary_files | wc"
        ls -lh gene_summary_files | wc

        # Run Python splice junction discovery script
        echo "python /process_coverage_info.py --label ${label} --cutoff ${cutoff} --output /cromwell_root/ --gene gene_summary_files"
        python /process_coverage_info.py --label ${label} --cutoff ${cutoff} --output /cromwell_root/ --gene gene_summary_files --threads ${threads}

        echo "Output files in /cromwell_root/:"
        ls -lh /cromwell_root/

    >>>

    output {
		File per_sample_gene_fraction="${label}_gene_fraction_above_15_coverage_per_sample.tsv.gz"
    }

    runtime {
      	docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
        memory: "${memory}"
        disks: "local-disk " + diskSpace + " HDD"
        bootDiskSizeGb: "${bootSizeDiskGb}"
        preemptible: preemptible
    }
}

task ProcessCoverageInfoIntervalTask {
    File intervalSummaryFilesGzipped
    Int diskSpace
    String memory
    String label
    Int cutoff
    Int preemptible
    String bootSizeDiskGb
    String dockerVersion
    Int threads

    command <<<
    	# log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

        tar -xvzf ${intervalSummaryFilesGzipped}

        echo "ls -lh"
        ls -lh

        echo "mkdir output"
        mkdir output

        echo "ls -lh interval_summary_files | wc"
        ls -lh interval_summary_files | wc

        # Run Python splice junction discovery script
        echo "python /process_coverage_info.py --label ${label} --cutoff ${cutoff} --output /cromwell_root/ --interval interval_summary_files"
        python /process_coverage_info.py --label ${label} --cutoff ${cutoff} --output /cromwell_root/ --interval interval_summary_files --threads ${threads}

        echo "Output files in /cromwell_root/:"
        ls -lh /cromwell_root/

    >>>

    output {
        File per_sample_interval_fraction="${label}_interval_fraction_above_15_coverage_per_sample.tsv.gz"
    }

    runtime {
      	docker: "vanallenlab/depth_of_coverage_summary:${dockerVersion}"
        memory: "${memory}"
        disks: "local-disk " + diskSpace + " HDD"
        bootDiskSizeGb: "${bootSizeDiskGb}"
        preemptible: preemptible
    }
}