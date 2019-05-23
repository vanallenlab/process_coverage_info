# docker build -t "vanallenlab/depth_of_coverage_summary:{}.0"

FROM vanallenlab/miniconda:3.6

COPY process_coverage_info.py /
COPY differential_coverage_analysis.py /

RUN pip install matplotlib
RUN pip install pandas
RUN pip install numpy
RUN pip install argparse
RUN pip install scipy

RUN apt-get clean
RUN apt-get install tabix
RUN ulimit -Hn 10240