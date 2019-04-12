FROM vanallenlab/miniconda:3.6

COPY process_coverage_info.py /

RUN pip install matplotlib
RUN pip install pandas
RUN pip install numpy
RUN pip install argparse