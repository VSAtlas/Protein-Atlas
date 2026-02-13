#!/usr/bin/env python
import argparse
import json
import os
import pickle
import random
import multiprocessing as mp

import numpy as np

from DeepCoy_duds.DeepCoy import DenseGGNNChemModel
from DeepCoy_duds.data_augmentation import construct_incremental_graph_freqs
from DeepCoy_duds import utils

def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess DeepCoy JSON datasets into cached pickle files."
    )
    parser.add_argument("--in-json", required=True, help="Input JSON dataset file.")
    parser.add_argument("--out-cache", required=True, help="Output cache pickle file.")
    parser.add_argument(
        "--dataset", default="zinc", help="Dataset name (default: zinc)."
    )
    parser.add_argument(
        "--is-training", action="store_true", help="Mark data as training data."
    )
    parser.add_argument(
        "--restrict", type=int, default=None, help="Only preprocess first N molecules."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes (default: 1).",
    )
    return parser.parse_args()


def compute_meta(raw_data, params):
    num_fwd_edge_types = len(utils.bond_dict) - 1
    num_edge_types = num_fwd_edge_types * (1 if params["tie_fwd_bkwd"] else 2)
    max_num_vertices = 0
    for g in raw_data:
        max_in = max([v for e in g["graph_in"] for v in [e[0], e[2]]], default=0)
        max_out = max([v for e in g["graph_out"] for v in [e[0], e[2]]], default=0)
        max_num_vertices = max(max_num_vertices, max_in, max_out)
    annotation_size = len(raw_data[0]["node_features_in"][0])
    return max_num_vertices, num_edge_types, annotation_size


def load_freq_dict(params):
    freq_dict = {}
    freq_file = params.get("subgraph_freq_file")
    if freq_file:
        with open(freq_file, "rb") as f:
            freq_dict = pickle.load(f)
    return freq_dict


def _process_single_entry(args):
    d, bucket_sizes, params, is_training_data, freq_dict = args
    bucket_sizes = np.array(bucket_sizes)
    incremental_results = []
    new_raw_data = []

    out_direc = "out"
    if not params["path_random_order"]:
        if params["multi_bfs_path"]:
            list_of_starting_idx = list(range(params["bfs_path_count"]))
        else:
            list_of_starting_idx = [0]
    else:
        node_length = len(d["node_features_" + out_direc])
        if params["multi_bfs_path"]:
            list_of_starting_idx = np.random.choice(
                node_length, params["bfs_path_count"], replace=True
            )
        else:
            list_of_starting_idx = [random.choice(list(range(node_length)))]

    for list_idx, starting_idx in enumerate(list_of_starting_idx):
        chosen_bucket_idx = np.argmax(
            bucket_sizes
            > max(
                max([v for e in d["graph_out"] for v in [e[0], e[2]]]),
                max([v for e in d["graph_in"] for v in [e[0], e[2]]]),
            )
        )
        chosen_bucket_size = bucket_sizes[chosen_bucket_idx]

        nodes_no_master = d["node_features_" + out_direc]
        edges_no_master = d["graph_" + out_direc]

        (
            incremental_adj_mat,
            distance_to_others,
            node_sequence,
            edge_type_masks,
            edge_type_labels,
            local_stop,
            edge_masks,
            edge_labels,
            overlapped_edge_features,
            new_compound_frequencies,
            new_compound_frequencies_edge,
        ) = construct_incremental_graph_freqs(
            params["dataset"],
            edges_no_master,
            chosen_bucket_size,
            len(nodes_no_master),
            nodes_no_master,
            params,
            is_training_data,
            freq_dict,
            initial_idx=starting_idx,
        )

        if params["sample_transition"] and list_idx > 0:
            incremental_results[-1] = [
                x + y
                for x, y in zip(
                    incremental_results[-1],
                    [
                        incremental_adj_mat,
                        distance_to_others,
                        node_sequence,
                        edge_type_masks,
                        edge_type_labels,
                        local_stop,
                        edge_masks,
                        edge_labels,
                        overlapped_edge_features,
                        new_compound_frequencies,
                        new_compound_frequencies_edge,
                    ],
                )
            ]
        else:
            incremental_results.append(
                [
                    incremental_adj_mat,
                    distance_to_others,
                    node_sequence,
                    edge_type_masks,
                    edge_type_labels,
                    local_stop,
                    edge_masks,
                    edge_labels,
                    overlapped_edge_features,
                    new_compound_frequencies,
                    new_compound_frequencies_edge,
                ]
            )
            new_raw_data.append(d)

    return incremental_results, new_raw_data


class PreprocessModel(DenseGGNNChemModel):
    def __init__(self, params, num_edge_types, freq_dict, workers):
        self.params = params
        self.num_edge_types = num_edge_types
        self.freq_dict = freq_dict
        self.workers = workers

    def calculate_incremental_results(
        self, raw_data, bucket_sizes, file_name, is_training_data
    ):
        if self.workers <= 1:
            return super().calculate_incremental_results(
                raw_data, bucket_sizes, file_name, is_training_data
            )

        if self.params.get("path_random_order"):
            print(
                "Warning: path_random_order enabled; falling back to single-worker preprocessing."
            )
            return super().calculate_incremental_results(
                raw_data, bucket_sizes, file_name, is_training_data
            )

        incremental_results = [[], []]
        new_raw_data = []
        tasks = [
            (d, bucket_sizes, self.params, is_training_data, self.freq_dict)
            for d in raw_data
        ]
        with mp.Pool(processes=self.workers) as pool:
            for idx, (incremental_list, raw_list) in enumerate(
                pool.imap(_process_single_entry, tasks)
            ):
                incremental_results[1].extend(incremental_list)
                new_raw_data.extend(raw_list)
                if idx % 50 == 0:
                    print(
                        "Finished calculating %d incremental matrices" % idx, end="\r"
                    )
        return incremental_results, new_raw_data


def main():
    args = parse_args()

    in_json = os.path.abspath(args.in_json)
    out_cache = os.path.abspath(args.out_cache)

    if not os.path.isfile(in_json):
        raise FileNotFoundError("Input JSON not found: %s" % in_json)

    with open(in_json, "r") as f:
        raw_data = json.load(f)

    if args.restrict is not None and args.restrict > 0:
        raw_data = raw_data[: args.restrict]

    if not raw_data:
        raise ValueError("No data loaded from %s" % in_json)

    params = DenseGGNNChemModel.default_params()
    params["dataset"] = args.dataset

    freq_dict = load_freq_dict(params)
    max_num_vertices, num_edge_types, annotation_size = compute_meta(raw_data, params)

    model = PreprocessModel(params, num_edge_types, freq_dict, args.workers)
    processed = model.process_raw_graphs(raw_data, args.is_training, in_json)

    payload = {
        "file_name": in_json,
        "dataset": params["dataset"],
        "is_training_data": bool(args.is_training),
        "processed": processed,
        "meta": {
            "max_num_vertices": max_num_vertices,
            "num_edge_types": num_edge_types,
            "annotation_size": annotation_size,
        },
    }

    out_dir = os.path.dirname(out_cache)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_cache, "wb") as f:
        pickle.dump(payload, f, pickle.HIGHEST_PROTOCOL)

    print("Wrote preprocessed cache to %s" % out_cache)


if __name__ == "__main__":
    main()
