from collections import namedtuple
import itertools
from typing import Dict, List

from dask.distributed import Future, get_client, rejoin, secede, as_completed, wait
import pandas as pd


ParamArg = namedtuple("ParamArg", ["inputs", "params"])


def task_factory(
    task_schema,
    client,
):
    task_gathered = {}
    for task in task_schema:
        
        # task["client_method"] = "submit"
        task_futures = getattr(client, task["client_method"])(
            task["callable"],
            *parse_inputs(task["inputs"],
                          task_gathered,
                          pd.json_normalize(task_schema)),
        )
        
        task_gathered[task["taskID"]] = task_futures

    return task_gathered[task["taskID"]]

def parse_inputs(
    inputs,
    task_futures: Dict[str, Future],
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
    elif isinstance(inputs, Future):
        return inputs
    else:
        # warnings.warn(f"Unknown input type: {type(inputs)}. "
        #               "Please use schema field to specify input types.")
        return inputs

def sub_task(arg_client=None):
    # Not advised to use this method
    client = arg_client if arg_client is not None else get_client()
    _seceded = False
    if task["gather"] is True:
        secede()
        task_gathered[task["taskID"]] = client.gather(
            task_futures
        )
        rejoin()
    else:
        task_gathered[task["taskID"]] = task_futures