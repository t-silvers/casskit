import itertools
from typing import Dict

from dask.distributed import Future, get_client, rejoin, secede
import networkx as nx
import numpy as np
import pandas as pd


class AssociationTasks:
    def __init__(
        self,
        name,
        description,
        inputs,
        tasks,
    ):
        self.name = name
        self.description = description
        self.inputs = inputs
        self.tasks = tasks

    @property
    def schema(self):
        return {
            "name": self.name,
            "description": self.description,
            "inputs": self.inputs,
            "tasks": self.tasks,
        }

    def make_tasks(self):
        pass

    def execute(self, inputs):
        pass

    def __repr__(self):
        return f"AssociationTasks(name={self.name}, description={self.description})"

    def task_factory(
        task_id,
        model_frame,
        design_response_ixs,
        space,
        ignore_support,
        groups_idxs,
    ):
        # Set up the computing environment
        client = get_client()

        # Prepare data
        _res0_future = client.submit(_task0, model_frame, design_response_ixs)
        _res1_future = client.submit(_task1, _res0_future)

        # Optimize
        _res2_futures = client.map(
            _task2, space, list(itertools.repeat(_res0_future, len(space))),
            list(itertools.repeat(groups_idxs, len(space)))
        )
        
        secede()

        _res3_futures = client.map(
            _task3,
            list(itertools.repeat(_res0_future, len(space))),
            _res2_futures,
            list(itertools.repeat(_res1_future, len(space))),
            list(itertools.repeat(ignore_support, len(space))),
        )
        
        _res3_results = client.gather(_res3_futures)

        rejoin()

        # Complete
        _res4_future = client.submit(_task4, space, _res3_results)
        _opt_fit = client.submit(_task2, _res4_future, _res0_future, groups_idxs)

        secede()
        opt_param = client.gather(_opt_fit)
        
        rejoin()
        
        return task_id, opt_param


# =================================================


def sample_log_liks(X, y, estimator, scale):
    y_pred = estimator.predict_expected(X) # for linear GLM, just dot product + intercept
    # y_pred = estimator.predict(X)
    return glm_log_liks.gaussian(y_pred=y_pred, y_true=y, scale=scale)

def calc_neg_ebic(X, y, estimator, scale, ignore_support=0):
    """Calculate the extended BIC for a given estimator."""
    log_lik = sample_log_liks(X, y, estimator, scale).sum()
    
    n_support = np.count_nonzero(estimator.coef_) - ignore_support
    
    return -info_criteria.ebic(log_lik=log_lik,
                               n_samples=X.shape[0],
                               n_features=X.shape[1],
                               n_support=n_support,
                               gamma='default',
                               fit_intercept=estimator.fit_intercept)

def _feature_selector(model_frame, ixs):
    response_ix, design_ix = ixs
    X = model_frame.iloc[:, design_ix].fillna(0)
    y = model_frame.iloc[:, response_ix]
    return X.values, y.values

def _linreg_noise_var(_task0):
    return lin_reg_noise_var.ViaRidge().fit(*_task0).scale_

def _cneqtl_model(alpha, Xy, groups_idxs):
    
    pens = {'no_pen': range(0, 3),
            'sparse': range(3, Xy[0].shape[1])}

    # groups_idxs = mod_kwargs.get("groups_idxs")

    mod = Glm(loss="lin_reg",
              penalty=SeparableSum(groups=pens,
                                   no_pen=NoPenalty(),
                                   sparse=SparseGroupLasso(pen_val=alpha,
                                                           groups=groups_idxs)),
              fit_intercept=False,
              standardize=False)

    return mod.fit(*Xy)

def _score_ebic(_task0, _task2_out, _task1_out, ignore_support):
    return calc_neg_ebic(*_task0, _task2_out, _task1_out, ignore_support)

def _eval_ebic(space, _task3_out):
    return space[np.argmax(_task3_out)]



# =================================================

from collections import namedtuple

model_frame = pd.DataFrame()
design_response_ixs = (0, 1)
groups_idxs = [range(1)]
space = [0.1, 0.2, 0.3]
ignore_support = 1

ParamArg = namedtuple("ParamArg", ["inputs", "params"])


task_schema = [
    {
        "taskID": "taskID1",
        "callable": _feature_selector,
        "inputs": [model_frame, design_response_ixs],
        "params": {},
        "gather": False,
    },
    {
        "taskID": "taskID2",
        "callable": _linreg_noise_var,
        "inputs": ["taskID1"],
        "params": {},
        "gather": False,
    },    
    {
        "taskID": "taskID3",
        "callable": _cneqtl_model,
        "inputs": [space, ParamArg("taskID1", space), ParamArg(groups_idxs, space)],
        "params": {"space": space},
        "gather": False,
    },
    {
        "taskID": "taskID4",
        "callable": _score_ebic,
        "inputs": ["taskID1", "taskID3", "taskID2", ignore_support],
        "params": {},
        "gather": True,
    },
    {
        "taskID": "taskID5",
        "callable": _eval_ebic,
        "inputs": [space, "taskID4"],
        "params": {},
        "gather": False,
    },
    {
        "taskID": "taskID6",
        "callable": _cneqtl_model,
        "inputs": ["taskID5", "taskID1", groups_idxs],
        "params": {},
        "gather": True,
    },
]


schema = pd.json_normalize(task_schema)

def parse_tasks(task, task_futures, schema):
    if len(task["params"]) == 0:
        return client.submit(
            task["callable"],
            *parse_inputs(task["inputs"], task_futures, schema),
        )
    
    else:
        return client.map(
            task["callable"],
            *parse_inputs(task["inputs"], task_futures, schema)
        )

def parse_inputs(
    inputs,
    tasks: Dict[Future],
    schema: pd.DataFrame
):
    if isinstance(inputs, list):
        return [parse_inputs(x) for x in inputs]
    elif isinstance(inputs, str):
        if inputs.startswith("taskID"):
            return tasks[schema.query("taskID == @inputs").iloc[0].taskID]
        else:
            return
    elif isinstance(inputs, ParamArg):
        return list(itertools.repeat(
            parse_inputs(getattr(inputs, "inputs")),
            len(getattr(inputs, "params"))
        ))


client = get_client()

task_futures = {}
task_gathered = {}

_seceded = False

for task in task_schema:

    if _seceded is True:
        if len(task["params"]) == 0:
            rejoin()
            _seceded = False
    
    task_futures = parse_tasks(
        task, task_gathered, schema
    )
    
    if task["gather"] is True:
        task_gathered[task["taskID"]] = client.gather(
            task_futures
        )
    else:
        task_gathered[task["taskID"]] = task_futures
    
    if len(task["params"]) > 0:
        if _seceded is False:
            secede()
            _seceded = True





G = nx.DiGraph()