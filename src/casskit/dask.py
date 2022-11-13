from contextlib import contextmanager
import os
from pathlib import Path
import subprocess

from dask_jobqueue import SLURMCluster
from dask.distributed import Client


class DaskCluster:
    def __init__(self, cores, memory, numworkers, time_limit):
        self.cores = cores
        self.memory = memory
        self.numworkers = numworkers
        self.time_limit = time_limit
        if time_limit is None:
            self.time_limit = self.timelimit()

    @classmethod
    @contextmanager
    def dask_cluster(cls, cores, memory, numworkers, time_limit):
        """Context manager to launch a Dask cluster on SLURM.

        Example:
        -------
        with DaskCluster(4, '4GB', 2, "0:30:00") as (cluster, client):
            # do stuff
        """
        cluster, client = cls(cores, memory, numworkers, time_limit).acquire_cluster()
        try:
            yield cluster, client

        finally:
            client.close()
            cluster.close()

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
            processes=1,
            queue='hbfraser,hns',
            memory=self.memory,
            walltime=self.time_limit,
            job_extra_directives=[
                '--job-name="smk-dask-worker"',
                '--propagate=NONE'
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
        return subprocess.check_output("TIME=$(squeue -j $SLURM_JOB_ID -h --Format TimeLimit); echo -n $TIME", shell=True, text=True)