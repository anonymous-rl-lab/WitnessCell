FROM python:3.12.14-slim
WORKDIR /opt/witnesscell
COPY requirements*.txt ./
RUN python -m pip install --no-cache-dir -r requirements.txt
COPY . .
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
ENTRYPOINT ["python", "-B", "reproduce.py"]
CMD ["--out", "/results", "--jobs", "2", "--compatible"]
