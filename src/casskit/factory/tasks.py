from collections import namedtuple
import itertools
from typing import Dict, List
import warnings

import dask.distributed
from dask.distributed import get_client, rejoin, secede
import pandas as pd


ParamArg = namedtuple("ParamArg", ["inputs", "params"])


def task_factory(
    task_schema,
    client: dask.distributed.Client = None,
    dry_run: bool = False,
):
    if dry_run is True:
        dry_run_results = {}
        for task in task_schema:
            try:
                inputs = parse_inputs(task["inputs"],
                                      dry_run_results,
                                      pd.json_normalize(task_schema))
                
                dry_run_results[task["taskID"]] = task.get("callable")(*inputs)
            
            except Exception as e:
                print(f"Error in task {task['taskID']}: {e}")
                print(f"Task: {task}")
                return inputs
        
        return dry_run_results
        
    task_gathered = {}
    for task in task_schema:
        task_futures = getattr(client, task["client_method"])(
            task["callable"],
            *parse_inputs(task["inputs"],
                          task_gathered,
                          pd.json_normalize(task_schema)),
            pure=task.get("pure", True),
        )
        
        task_gathered[task["taskID"]] = task_futures

    return task_gathered[task["taskID"]]

def parse_inputs(
    inputs,
    task_futures: Dict[str, dask.distributed.Future],
    schema: pd.DataFrame
):
    if isinstance(inputs, list):
        return [parse_inputs(x, task_futures, schema) for x in inputs]
    elif isinstance(inputs, str):
        if inputs.startswith("taskID"):
            return task_futures[schema.query("taskID == @inputs").iloc[0].taskID]
        else:
            return inputs
    elif isinstance(inputs, ParamArg):
        return list(itertools.repeat(
            parse_inputs(getattr(inputs, "inputs"), task_futures, schema),
            len(getattr(inputs, "params"))
        ))
    elif isinstance(inputs, dask.distributed.Future):
        return inputs
    else:
        warnings.warn(f"Unknown input type: {type(inputs)}. "
                      "Please use schema field to specify input types.")
        return inputs

def sub_task(
    task: Dict,
    task_futures: Dict[str, dask.distributed.Future],
    arg_client: dask.distributed.Client = None,
):
    # Not advised to use this method
    client = arg_client if arg_client is not None else get_client()
    _seceded = False
    if task["gather"] is True:
        secede()
        task_futures[task["taskID"]] = client.gather(
            task_futures
        )
        rejoin()
    else:
        task_futures[task["taskID"]] = task_futures
