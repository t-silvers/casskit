from contextlib import contextmanager
import os
from pathlib import Path
import subprocess

from dask_jobqueue import SLURMCluster
from dask.distributed import Client

import casskit.config as config


class DaskCluster:
    def __init__(
        self,
        cores: int = 4,
        memory: str = "8GB",
        numworkers: int = 5,
        threads: int = 4,
        time_limit: str = None,
        log_dir: str = None,
    ):
        self.cores = cores
        self.memory = memory
        self.numworkers = numworkers
        self.threads = threads
        self.time_limit = time_limit
        
        if time_limit is None:
            self.time_limit = self.timelimit()
        
        elif time_limit == "default":
            self.time_limit = "0:30:00"
        
        self.log_dir = log_dir
        if log_dir is None:
            self.log_dir = config.get_cache() / "dask_logs"
            self.log_dir.mkdir(exist_ok=True)

    @classmethod
    @contextmanager
    def dask_cluster(
        cls,
        cores: int = 4,
        memory: str = "8GB",
        numworkers: int = 5,
        threads: int = 4,
        time_limit: str = None
    ):
        """Context manager to launch a Dask cluster on SLURM.

        Example:
        -------
        with DaskCluster(4, '4GB', 2, "0:30:00") as (cluster, client):
            # do stuff
        """
        cluster, client = cls.start(cores, memory, numworkers, threads, time_limit)
        try:
            yield cluster, client

        finally:
            client.close()
            cluster.close()

    @classmethod
    def start(
        cls,
        cores: int = 4,
        memory: str = "8GB",
        numworkers: int = 5,
        threads: int = 4,
        time_limit: str = None
    ):
        return cls(cores, memory, numworkers, threads, time_limit).acquire_cluster()

    def acquire_cluster(self):
        """
        see here https://github.com/dask/distributed/issues/2694
        dask.config.set({'distributed.scheduler.allowed-failures': 10})

        Check nanny's memory allocation threshold:
        str(dask.config.get("distributed.nanny.environ.MALLOC_TRIM_THRESHOLD_")) # '65536'
        """    

        os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'
        
        cluster = SLURMCluster(
            cores=self.cores,
            processes=self.threads,
            queue='hbfraser,hns',
            memory=self.memory,
            walltime=self.time_limit,
            job_extra_directives=[
                '--job-name="smk-dask-worker"',
                '--propagate=NONE',
                f'--error={self.log_dir}/smk-dask-worker-%j.err',
                f'--output={self.log_dir}/smk-dask-worker-%j.out'
            ],
            job_script_prologue=[
                'export MALLOC_TRIM_THRESHOLD_=0',
                'export PYTHONPATH=$PYTHONPATH:/workflow'
            ],
            worker_extra_args=[
                '--no-dashboard'
            ]
        )

        cluster.scale(self.numworkers)
        client = Client(cluster, timeout="30s")
        client.wait_for_workers(n_workers=self.numworkers)
        
        return cluster, client

    @staticmethod
    def timelimit():
        SLURM_JOB_TIME_LIMIT = (
            "TIME=$(squeue -j $SLURM_JOB_ID -h --Format TimeLimit); echo -n $TIME"
        )
        return subprocess.check_output(SLURM_JOB_TIME_LIMIT,
                                       shell=True,
                                       text=True)