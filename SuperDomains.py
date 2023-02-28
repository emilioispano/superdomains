#!/usr/bin/python3

from enum import Enum
from multiprocessing import Manager, Process
from rdflib import Graph, Namespace, OWL
from math import log
import io
from ctypes import CDLL, c_int, POINTER
import numpy as np
from typing import Set
import math
from decimal import Decimal
import argparse
import logging
from py4j import java_gateway
import pickle
from concurrent.futures import ThreadPoolExecutor
import random
import subprocess
import networkx as nx
from networkx.readwrite import gml
from pymongo import MongoClient
import multiprocessing
import os
import sys
import shutil

global iternum

gateway = java_gateway.JavaGateway().launch_gateway(port=0, classpath='jgrapht-core-1.3.0.jar')
ConnectivityInspector = gateway.jvm.org.jgrapht.alg.connettivity.ConnectivityInspector
GraphTypeBuilder = gateway.jvm.org.jgrapht.graph.builder.GraphTypeBuilder
DefaultWeightedEdge = gateway.jvm.org.jgrapht.graph.DefaultWeightedEdge
SimpleWeightedGraph = gateway.jvm.org.jgrapht.graph.SimpleWeightedGraph
AsSubgraph = gateway.jvm.org.jgrapht.graph.AsSubgraph
AsWeightedGraph = gateway.jvm.org.jgrapht.graph.AsWeightedGraph
SpanningTreeAlgorithm = gateway.jvm.org.jgrapht.alg.interfaces.SpanningTreeAlgorithm


def read_dom_file(file):
    """
    Reads a file containing domain information and returns a dictionary mapping a protein ID to a dictionary containing
    information about its domains.

    :param file: the file to read
    :return: a dictionary mapping a protein ID to a dictionary containing information about its domains
    """
    mapp = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            data = line.split("\t")
            if len(data) < 12:
                continue
            pid = data[0]
            ipr = data[11]
            p = int(data[6])
            if pid not in mapp:
                mapp[pid] = {}
            if p in mapp[pid]:
                mapp[pid][p].add(ipr)
            else:
                mapp[pid][p] = {ipr}
    return mapp


def parse_blast_file(file, setting):
    """
    Parses a BLAST file and returns a dictionary mapping a template sequence ID to a set of subject sequence IDs.

    :param file: the BLAST file to parse
    :param setting: the settings to use
    :return: a dictionary mapping a template sequence ID to a set of subject sequence IDs
    """
    blast_map = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("T"):
                continue
            data = line.split()
            template = data[0]
            e_value = -log(float(data[10]))
            if e_value >= setting.blast_thr:
                if template not in blast_map:
                    blast_map[template] = set()
                blast_map[template].add(data[1])
    return blast_map


def merge_thread_results(domains_thread_results):
    results = {}

    for dt in domains_thread_results:
        try:
            results.update(dt.output)
        except FileNotFoundError as ex:
            logging.getLogger('Utilities').error(ex)

    return results


def delete_directory(path):
    """
    Deletes a directory and all its contents.

    :param path: the path of the directory to delete
    """
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(path)


def clean_directory(path, start_dir):
    """
    Deletes all the contents of a directory except the starting directory itself.
    """
    for root, dirs, files in os.walk(path):
        for file in files:
            os.remove(os.path.join(root, file))
        for dirr in dirs:
            if os.path.abspath(os.path.join(root, dirr)) != start_dir:
                shutil.rmtree(os.path.join(root, dirr))


def delete_file(file):
    os.remove(file)


def create_directory(dirr):
    os.makedirs(dirr, exist_ok=True)


def create_file(file):
    open(file, 'a').close()


def directory_exists(dirr):
    return os.path.exists(dirr)


def file_exists(file):
    return os.path.exists(file)


def sort_uids(mapp):
    tmp = {}
    uids = []

    for uid, t in mapp.items():
        if t in tmp:
            tmp[t].append(uid)
        else:
            tmp[t] = [uid]

    keys = list(tmp.keys())
    keys.sort(reverse=True)

    for v in keys:
        for t in tmp[v]:
            uids.append(t)

    return uids


def read_clustering_data(file, database, setting):
    clusters = {}
    prot_ids = set()
    go_ids = set()
    goa_coll = database.get_collection(setting['collectionGOA'])

    if os.path.exists(file):
        with open(file) as f:
            next(f)
            repr_prot = ""
            it = 0

            for line in f:
                data = line.strip().split()

                if it > 0 and repr_prot != data[0]:
                    if repr_prot not in clusters:
                        clusters[repr_prot] = {'repr': repr_prot, 'protIds': prot_ids.copy(), 'goIds': go_ids.copy()}
                        prot_ids.clear()
                        go_ids.clear()

                repr_prot = data[0]
                prot_ids.add(data[1])

                for doc in goa_coll.find({'uid': data[1]}):
                    go_ids.add(doc['goid'])

                it += 1

    return clusters


def get_repr_to_removes(repr_score, setting):
    if len(repr_score.keys()) > setting['maxNumRepr']:
        sorted_repr = sort_uids(repr_score)
        return set(sorted_repr[:len(sorted_repr) - setting['maxNumRepr']])
    else:
        return set()


def convert_sec_to_day(n):
    start_value = n

    days = n // (24 * 3600)

    n %= (24 * 3600)
    hours = n // 3600

    n %= 3600
    minutes = n // 60

    n %= 60
    seconds = n

    print(f"{days} d {hours} h {minutes} m {seconds} s ({start_value} s)")


class Node:
    def __init__(self):
        self.ontid = ''
        self.parents = set()
        self.children = set()
        self.parent_nodes = {}
        self.child_nodes = {}
        self.id = 0
        self.flag = False
        self.defn = None
        self.comment = None
        self.namespace = None
        self.disjoint = set()
        self.int_disjoint = False
        self.namespace = None
        self.weight = 0.0
        self.abs_value = 0
        self.freq = -1
        self.real_freq = 0.0
        self.real_count = 0
        self.int_pid = set()
        self.z_score = -10

    def set_weight(self, weight):
        self.weight = weight

    def get_weight(self):
        return self.weight

    def get_round_weight(self, round_to):
        if round_to < 0:
            round_to = 0

        bd = Decimal(self.weight)
        bd = bd.quantize(Decimal('1.' + '0' * round_to), rounding='ROUND_HALF_UP')
        w = bd.to_eng_string()
        return float(w)

    def add_weight(self, w):
        self.weight += w

    def get_round_z_score(self, round_to):
        if round_to < 0:
            round_to = 0

        bd = Decimal(self.z_score)
        bd = bd.quantize(Decimal('1.' + '0' * round_to), rounding='ROUND_HALF_UP')
        z = bd.to_eng_string()
        return float(z)

    def add_abs_value(self, n):
        self.abs_value += n

    def get_abs_value(self):
        return self.abs_value

    def set_frequency(self, freq):
        self.freq = freq

    def get_frequency(self):
        return self.freq

    def get_ic(self):
        if self.abs_value == 0:
            return 0
        return -math.log(self.freq)

    def set_real_frequency(self, freq):
        self.real_freq = freq

    def get_real_frequency(self):
        return self.real_freq

    def set_real_count(self, real_count):
        self.real_count = real_count

    def get_real_count(self):
        return self.real_count

    def add_protein_id(self, prot_id):
        self.int_pid.add(prot_id)

    def get_protein_id_list(self):
        return self.int_pid

    def has_gos(self):
        return len(self.int_pid) > 0

    def set_ont_id(self, ontid):
        self.ontid = ontid

    def get_ont_id(self):
        return self.ontid

    def equals(self, node):
        return self.ontid == node.get_ont_id()

    def add_parent(self, edge):
        self.parents.add(edge)

    def add_child(self, edge):
        self.children.add(edge)

    def add_parent_node(self, n, e):
        self.parent_nodes[n] = e

    def add_child_node(self, n, e):
        self.child_nodes[n] = e

    def get_parents(self):
        return self.parents

    def get_children(self):
        return self.children

    def contain_parent(self, n):
        return n in self.parent_nodes

    def contain_child(self, n):
        return n in self.child_nodes

    def is_leaf(self):
        return len(self.children) == 0

    def set_flag(self, flag):
        self.flag = flag

    def is_flagged(self):
        return self.flag

    def remove_child(self, e):
        self.children.remove(e)

    def remove_parent(self, e):
        self.parents.remove(e)

    def remove_parent_node(self, n):
        p = self.parent_nodes[n]
        self.remove_parent(p)
        self.parent_nodes.pop(n)

    def remove_child_node(self, n):
        c = self.child_nodes[n]
        self.remove_child(c)
        self.child_nodes.pop(n)

    def is_root(self):
        return len(self.parents) == 0

    def clean(self):
        self.flag = False
        self.weight = 0.0
        self.z_score = -10

    def get_id(self):
        return self.id

    def set_id(self, idd):
        self.id = idd

    def set_namespace(self, namespace):
        if namespace == "biological_process":
            self.namespace = OntNamespace.biological_process
        elif namespace == "molecular_function":
            self.namespace = OntNamespace.molecular_function
        elif namespace == "cellular_component":
            self.namespace = OntNamespace.cellular_component

    def get_namespace(self):
        return self.namespace

    def set_definition(self, defn):
        self.defn = defn

    def set_comment(self, comment):
        self.comment = comment

    def get_definition(self):
        return self.defn

    def get_comment(self):
        return self.comment

    def add_disjoint(self, goid):
        if isinstance(goid, set):
            self.disjoint.update(goid)
        else:
            self.disjoint.add(goid)

    def has_disjoint(self):
        return bool(self.disjoint)

    def is_disjoint_with(self, goid):
        return goid in self.disjoint

    def set_int_disjoint(self, disj):
        self.int_disjoint = disj

    def is_int_disjoint(self):
        return self.int_disjoint

    @staticmethod
    def is_go_node(node):
        return isinstance(node, Node)


class AlignmentResults:
    def __init__(self):
        self.map = {}

    def get_map(self):
        return self.map

    def add_result(self, uid, sim):
        self.map[uid] = sim

    def get_uid_over(self, thr):
        tmp = {}
        for uid, sim in self.map.items():
            if sim >= thr:
                tmp[uid] = sim
        return tmp

    def get_uid_over_concurrent(self, thr):
        tmp = {}
        with ThreadPoolExecutor() as executor:
            for uid, sim in executor.map(lambda x: (x[0], x[1]), self.map.items()):
                if sim >= thr:
                    tmp[uid] = sim
        return tmp

    def get_uid(self):
        return self.map


class ArgotGraph:
    def __init__(self, uidscore, pid, database, blast_prot, thread_id):
        self.settings = Settings.init()
        self.filein = None
        self.retparam = None
        self.work_dir = self.settings.work_dir
        self.iternum = 0
        self.sg = None
        self.database = database
        self.buffer = set()
        self.randomseed = None
        self.cluster = {}
        self.prot_clusters = {}
        self.uids = []
        self.cluster = {}
        self.randomseed = self.get_random_seed()
        self.retparam = RetrieveProteins.init(True)
        self.repr_prot_ids = []
        self.repr_score = {}
        self.outfile = None

        self.work_dir = self.settings.get_work_dir() + "/" + str(thread_id)
        file_in_clust = self.work_dir + "/" + pid + "_" + self.randomseed + ".fa.tmp"
        self.filein = self.work_dir + "/" + pid + "_" + self.randomseed + ".fa"
        fileq = self.filein + "_" + str(thread_id) + ".q"
        filenwscore = self.filein + ".nwscore"
        tmp_clust = os.path.join(self.work_dir, "mmclust_tmp")
        out_clust = self.work_dir + "/" + pid + "_output"

        create_directory(self.work_dir)

        if blast_prot is not None and len(blast_prot) > 0:
            for prot in blast_prot:
                uidscore.pop(prot, None)

        if len(uidscore.keys()) > self.settings.get_max_seqs():
            uids = sort_uids(uidscore)
        else:
            uids = list(uidscore.keys())

        self.retparam.write_prot_file(uids, file_in_clust, pid, database)

        self.run_mmlin_clust(file_in_clust, out_clust, tmp_clust)
        self.prot_clusters = read_clustering_data(out_clust + "_cluster.tsv", self.database, self.settings)

        if len(self.prot_clusters.keys()) > 0:
            for reprId in self.prot_clusters.keys():
                self.repr_score[reprId] = uidscore[reprId]

            delete_file(file_in_clust)

            repr_rot_ids = list(self.prot_clusters.keys())
            self.retparam.write_prot_file_repr(repr_rot_ids, self.filein, pid, database)

            with open(filenwscore, "w") as outputnwscore:
                for p in uids:
                    outputnwscore.write(p + "\t" + str(uidscore[p]) + "\n")

            if not os.path.exists(tmp_clust):
                os.makedirs(tmp_clust)

            with open(self.settings.get_input_fasta(), "rb") as rafinput, open(fileq, "w") as output:
                self.retparam.write_query(rafinput, output, pid)

            if self.settings.is_rom_used():
                outfile = self.filein + "_" + str(thread_id) + ".aln"
                self.run_glsearch_alignment(fileq, self.filein, outfile)
                b8 = Blast8(outfile, self.settings.get_similarity_thr(), self.sg)
            else:
                lir = LoadInRam()
                self.run_glsearch_alignment_lir(fileq, self.filein, lir)
                b8 = Blast8(lir, self.settings.get_similarity_thr(), self.sg)
                lir.close()

            self.sg = b8.get_graph()
            self.targets = b8.get_targets()

            if self.settings.get_max_iteration() > 0 and len(self.targets) > 0:
                dirr = os.path.join(self.work_dir, str(self.randomseed))
                os.mkdir(dirr)
                self.deep_search(dirr, self.targets, self.database, self.settings)
                for uid in self.modularity(self.targets, self.sg):
                    score = 0.0
                    graphpath = nx.shortest_path(self.sg, source=uid, target=pid, weight='weight')
                    ledge = [self.sg[u][v]['weight'] for u, v in zip(graphpath, graphpath[1:])]
                    for edge in ledge:
                        score += self.to_score(edge)
                    score /= pow(2, len(ledge) - 1)
                    self.cluster[uid] = score
            else:
                if pid in self.sg:
                    for edge in self.sg.edges(pid):
                        self.cluster[self.sg[edge[1]]] = self.to_score((self.sg[edge]))

        if not self.settings.is_rom_used():
            os.remove(fileq)

            if self.outfile is not None:
                os.remove(self.outfile)

            os.remove(self.filein)

    def deep_search(self, dirr, targets, database, settings):
        global iternum
        deeptargets = set()
        iternum += 1

        lir = None
        if not settings.is_rom_used():
            lir = LoadInRam()

        for tid in targets:
            tidq = os.path.join(dirr, tid)
            self.retparam.write_single_prot(tid, tidq, database)
            tmp = os.path.join(self.work_dir, "mmseq_tmp")
            if not os.path.exists(tmp):
                os.mkdir(tmp)

            if settings.is_rom_used():
                outfile = os.path.join(dirr, tid + ".aln")
                self.run_mmsearch_alignment(tidq, self.filein, outfile, tmp)
                b8 = Blast8(outfile, settings.get_similarity_thr(), self.sg)
            else:
                self.run_mmsearch_alignment_lir(tidq, self.filein, lir, tmp)
                b8 = Blast8(lir, settings.get_similarity_thr(), self.sg)
                if lir:
                    lir.reset()

            self.sg = b8.get_graph()

            for uid in b8.get_targets():
                if uid not in self.buffer:
                    deeptargets.add(uid)

            self.buffer.update(deeptargets)

        if iternum <= settings.get_max_iteration() and deeptargets:
            if settings.isPrintgml():
                self.print_gml(iternum)

            self.deep_search(dirr, deeptargets, database, settings)

        if lir:
            lir.close()

        return self.buffer

    @staticmethod
    def modularity(pid, sg):
        ssn = SSNetwork(pid, sg)
        ids = ssn.get_cluster()
        tmp = set(ids)
        if pid in tmp:
            tmp.remove(pid)
        return tmp

    @staticmethod
    def to_score(evalue):
        if evalue == 0.0:
            return 300

        score = -math.log10(evalue)

        if score > 300:
            return 300

        return score

    def run_cmd(self, listargs, lir=None):
        pb = subprocess.Popen(listargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        if lir is None:
            logg = os.path.join(self.work_dir, "log")

            with open(logg, "wb") as f:
                for line in pb.stdout:
                    f.write(line)
        else:
            builder = []

            while True:
                line = pb.stdout.readline().decode()
                if not line:
                    break
                builder.append(line)

            lir.writeBytes(''.join(builder).encode())

        exitvalue = pb.wait()

        return exitvalue == 0

    def run_glsearch_alignment(self, fileq, filein, outfile, num_threads=1):
        selfalignargs = [
            self.settings.get_glsearch_alignment_cmd(),
            fileq,
            filein,
            "-m",
            "8",
            "-T",
            str(num_threads),
            "-O",
            outfile
        ]
        return self.run_cmd(selfalignargs)

    def run_glsearch_alignment_lir(self, fileq, filein, lir, num_threads=1):
        selfalignargs = [
            self.settings.get_glsearch_alignment_cmd(),
            "-m",
            "8",
            "-T",
            str(num_threads),
            fileq,
            filein
        ]
        return self.run_cmd(selfalignargs, lir)

    def run_mmsearch_alignment(self, fileq, filein, outfile, tmp_dir, num_threads=1):
        selfalignargs = [
            self.settings.get_mmsearch_alignment_cmd(),
            "easy-search",
            "--search-type",
            "1",
            "--threads",
            str(num_threads),
            fileq,
            filein,
            outfile,
            tmp_dir
        ]
        b = self.run_cmd(selfalignargs)
        if b:
            delete_directory(tmp_dir)
        return b

    def run_mmsearch_alignment_lir(self, fileq, filein, lir, tmp_dir, num_threads=1):
        selfalignargs = [
            self.settings.get_mmsearch_alignment_cmd(),
            "easy-search",
            "--search-type",
            "1",
            "--threads",
            str(num_threads),
            fileq,
            filein,
            tmp_dir
        ]
        b = self.run_cmd(selfalignargs, lir)
        if b:
            delete_directory(tmp_dir)
        return b

    def run_mmlin_clust(self, filein, outfile, tmp_dir, num_threads=1):
        selfalignargs = [
            self.settings.get_mmsearch_alignment_cmd(),
            "easy-linclust",
            "--min-seq-id",
            "0.90",
            "--threads",
            str(num_threads),
            filein,
            outfile,
            tmp_dir
        ]
        b = self.run_cmd(selfalignargs)
        if b:
            delete_directory(tmp_dir)
        return b

    @staticmethod
    def get_random_seed():
        saltchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
        salt = "".join(random.choice(saltchars) for _ in range(18))
        return salt

    def get_graph(self):
        return self.sg

    def print_gml(self, count):
        with open(f"{self.randomseed}_{count}.gml", "w") as out:
            gmlexp = gml.GmlExporter()
            gmlexp.add_graph(self.sg)
            out.write(gmlexp.dumps())

    def get_cluster(self):
        return self.cluster

    def get_prot_clusters(self):
        return self.prot_clusters


class Blast8:
    def __init__(self, file_or_lir, thr, sg):
        self.sg = sg
        self.targets = set()
        self.settings = Settings().init()
        self.counthit = 0
        self.bufferuid = set()

        if isinstance(file_or_lir, str):
            self.read_file(file_or_lir, thr)
        elif isinstance(file_or_lir, LoadInRam):
            self.read_data(file_or_lir, thr)
        else:
            raise ValueError("Expected either a file path or a LoadInRam instance")

    def read_file(self, file, thr):
        self.bufferuid = set()
        try:
            with open(file, 'r') as f:
                for line in f:
                    line = line.strip()
                    data = line.split('\t')
                    if data[0] != data[1]:
                        if not self.sg.__contains__(data[0]):
                            self.sg.add_vertex(data[0])

                        if not self.sg.__contains__(data[1]):
                            self.sg.add_vertex(data[1])

                        evalue = float(data[10])
                        score = -math.log(evalue)
                        if score >= thr:
                            if score >= self.settings.get_min_score():
                                self.counthit += 1
                            self.bufferuid.add(data[1])

                            if float(data[2]) <= self.settings.get_deep_thr():
                                self.targets.add(data[1])

                            if not self.sg.__contains_edge__(data[0], data[1]):
                                e = self.sg.add_edge(data[0], data[1])
                                self.sg.set_edge_weight(e, evalue)
                            else:
                                e = self.sg.get_edge(data[0], data[1])
                                if evalue < self.sg.get_edge_weight(e):
                                    self.sg.set_edge_weight(e, evalue)
        except FileNotFoundError as ex:
            logging.getLogger("Blast8").log(0, None, ex)
        except IOError as ie:
            logging.getLogger("Blast8").log(0, None, ie)

    def read_data(self, lir, thr):
        self.bufferuid = set()
        lines = lir.get_content().split('\n')

        for line in lines:
            data = line.split('\t')

            if len(data) > 1 and data[0] != data[1]:
                if not self.sg.__contains__(data[0]):
                    self.sg.add_vertex(data[0])

                if not self.sg.__contains__(data[1]):
                    self.sg.add_vertex(data[1])

                evalue = float(data[10])
                score = -math.log(evalue)

                if score >= thr:
                    if score >= self.settings.get_min_score():
                        self.counthit += 1
                    self.bufferuid.add(data[1])

                    if float(data[2]) <= self.settings.get_deep_thr():
                        self.targets.add(data[1])

                    if not self.sg.__contains_edge__(data[0], data[1]):
                        e = self.sg.add_edge(data[0], data[1])
                        self.sg.set_edge_weight(e, evalue)
                    else:
                        e = self.sg.get_edge(data[0], data[1])
                        if e is not None and evalue < self.sg.get_edge_weight(e):
                            self.sg.set_edge_weight(e, evalue)

    def get_graph(self):
        return self.sg

    def get_targets(self):
        return self.targets

    def to_go_deep(self):
        return self.counthit <= self.settings.get_hitnum()

    def get_uid(self) -> Set[str]:
        return self.bufferuid


class BMADistance:
    def __init__(self, distance_method):
        self.distance_method = distance_method

    def semantic_similarity(self, list1, list2, sim_limit=None):
        bcol = np.full(len(list2), self.distance_method.get_min_distance())
        brow = np.full(len(list1), self.distance_method.get_min_distance())

        for i, n in enumerate(list1):
            for j, m in enumerate(list2):
                dist = self.distance_method.compute_distance(n, m)
                brow[i] = max(brow[i], dist)
                bcol[j] = max(bcol[j], dist)

        sim_score = (np.sum(brow) + np.sum(bcol)) / (len(list1) + len(list2))

        if sim_limit is not None and sim_score < sim_limit:
            return 0.0
        else:
            return sim_score


class BMMDistance:
    def __init__(self, distance_method):
        self.distance_method = distance_method

    def semantic_similarity(self, list1, list2):
        simat = [[0.0 for _ in range(len(list2))] for _ in range(len(list1))]

        for i in range(len(list1)):
            n = list1[i]
            for j in range(len(list2)):
                simat[i][j] = self.distance_method.compute_distance(n, list2[j])

        k = 0
        tmp = 0.0
        totsim = 0.0
        m = min(len(list1), len(list2))
        for n in range(m):
            for i in range(k, len(list1)):
                tmp = max(simat[i][k], tmp)
            for j in range(k + 1, len(list2)):
                tmp = max(simat[k][j], tmp)
            k += 1
            totsim += tmp
            tmp = 0.0

        sim = totsim / max(len(list1), len(list2))
        return sim


class ClusterData:
    def __init__(self, representative, proteins, gos=None):
        self.representative = representative
        self.proteins = set(proteins)
        self.gos = set(gos) if gos else set()

    def get_representative(self):
        return self.representative

    def set_representative(self, representative):
        self.representative = representative

    def get_proteins(self):
        return self.proteins

    def set_proteins(self, proteins):
        self.proteins = set(proteins)

    def get_gos(self):
        return self.gos

    def set_gos(self, gos):
        self.gos = set(gos)

    def add_protein(self, prot_id):
        self.proteins.add(prot_id)

    def add_all_proteins(self, prot_ids):
        self.proteins.update(prot_ids)

    def add_go(self, go):
        self.gos.add(go)

    def add_all_gos(self, gos):
        self.gos.update(gos)


class DataDomains:
    def __init__(self, database, qdoms, uid=None):
        if not uid is None:
            self.database = database
            self.weightscore = {}
            self.aligners = AlignmentResults()
            tdoms = self.get_prot_doms(uid)
            self.weigth_score(tdoms)
            self.weigth_score(qdoms)
            nw = NW(qdoms, tdoms, self.weigthscore)
            self.aligners.add_result(uid, nw.get_similarity())
        else:
            self.database = database
            self.weigthscore = {}
            self.aligners = AlignmentResults()
            self.exit = False
            self.weigth_score(qdoms)
            seedoms = self.select_seed_dom(qdoms)
            inter = self.intersection(seedoms)

            for uid in inter:
                tdoms = self.get_prot_doms(uid)
                self.weigth_score(tdoms)
                sdom = True

                if len(seedoms) == 1:
                    q = qdoms[0]
                    for tdom in tdoms:
                        if q != tdom:
                            sdom = False
                            break

                if sdom and qdoms[0] == tdoms[0]:
                    self.aligners.add_result(uid, 1.0)
                else:
                    nw = NW(qdoms, tdoms, self.weigthscore)
                    self.aligners.add_result(uid, nw.get_similarity())

    def select_seed_dom(self, qdoms):
        numd = 1

        if len(qdoms) > 1:
            numd = len(qdoms) // 2

        seed = []
        mapp = {}
        domainfreq = []
        collection = self.database[Settings().init().get_collection_frequences()]

        for ipr in qdoms:
            freq = collection.find_one({"ipr": ipr})["freq"]
            mapp[freq] = ipr
            seed.append(freq)

        alg = Settings().init().get_seed_dom()

        if alg == "minfreq":
            seed.sort()
        elif alg == "maxfreq":
            seed.sort(reverse=True)
        elif alg == "all":
            numd = len(seed)

        for d in seed[:numd]:
            domainfreq.append(mapp[d])
        return domainfreq

    def intersection(self, doms):
        inter = set()
        list_ = []
        collection = self.database[Settings().init().get_collection_interpro()]

        for ipr in doms:
            buffer = set()
            for doc in collection.find({"ipr": ipr}):
                buffer.add(doc["uid"])
            list_.append(buffer)

        inter = list_[0]

        if len(list_) > 1:
            for i in range(1, len(list_)):
                inter &= list_[i]
        return inter

    def get_prot_doms(self, uid):
        mapp = {}
        collection = self.database.get_collection(Settings().init().get_collection_interpro())

        for doc in collection.find({"uid": uid}):
            ipr = doc["ipr"]
            pos_list = doc["pos"]

            for p in pos_list:
                if p in mapp:
                    mapp[p].add(ipr)
                else:
                    mapp[p] = {ipr}

        return [",".join(mapp[p]) for p in sorted(mapp.keys())]

    def weigth_score(self, doms):
        collection = self.database[Settings().init().get_collection_frequences()]

        for ipr in doms:
            if ipr not in self.weightscore:
                cursor = collection.find({"ipr": ipr})
                if cursor.count() > 0:
                    doc = cursor.next()
                    invfreq = 1 / doc["freq"]
                    weight = math.log(invfreq) / math.log(2)
                    self.weightscore[ipr] = weight
                else:
                    self.exit = True

    def get_exit(self):
        return self.exit

    def get_results(self):
        return self.aligners


class Direction(Enum):
    transitive = 1
    propagationUp = 2
    propagationDown = 3


class Disjoint:
    def __init__(self):
        self.graph = None
        self.relation = None
        self.visited = set()
        self.mapset = {}

    def add_disjoint(self, graph, bpic, ccic, mfic):
        self.relation = Relation().instance()
        self.visited = set()
        self.graph = graph

        try:
            print("Process")
            self.mapset = {}
            root = graph.get_gonode(graph.ProcessID())
            self.start(root, bpic)

            print("Function")
            root = graph.get_gonode(graph.FunctionID())
            self.start(root, ccic)

            print("Component")
            root = graph.get_gonode(graph.ComponentID())
            self.start(root, mfic)
        except NodeNotFoundException as ex:
            print('cosa ciera qui?')

    def start(self, root, ic):
        edges = self.graph.get_gochildren(root)
        for e in edges:
            if e.get_node().is_gonode():
                disjoint = self.find_disjoint(e.get_node(), ic)
                for goid, parents in disjoint.items():
                    self.graph.get_gonode(goid).add_disjoint(parents)

    def find_disjoint(self, node, ic):
        buffer = {}

        for e in self.graph.get_gochildren(node):
            if e.get_node().is_gonode():
                n = e.get_node()
                if 0 < n.get_ic() <= ic:
                    buffer[n.get_ontid()] = set()
                    self.walk_descendants(n, n.get_ontid(), buffer)

                    self.mapset.clear()
                    self.visited.clear()

        self.mapset.clear()
        self.visited.clear()

        return self.get_disjoint(buffer)

    def walk_descendants(self, n, parent, buffer):
        if n.get_ontid() in self.mapset:
            for goid in self.mapset[n.get_ontid()]:
                buffer[goid].add(parent)
            self.mapset[n.get_ontid()].add(parent)
        else:
            self.mapset[n.get_ontid()] = {parent}

        self.visited.add(n)

        if not n.is_leaf():
            for e in n.get_children():
                if isinstance(e.get_node(), Node):
                    if self.relation.get_relation(e.get_relationship_type()) != Direction.propagationUp:
                        self.walk_descendants(e.get_node(), parent, buffer)

    @staticmethod
    def get_disjoint(buffer):
        disjoint = {}

        markers = list(buffer.keys())
        for parent in markers:
            parents = buffer[parent]
            if parents:
                disjoint[parent] = set()
                for marker in markers:
                    if marker not in parents:
                        disjoint[parent].add(marker)
            else:
                disjoint[parent] = set(markers) - {parent}

        return disjoint

    @staticmethod
    def write_disjoint(disjoint):
        for goid, parents in disjoint.items():
            print(goid, end='')
            for marker in parents:
                print('\t' + marker, end='')
            print()


class DomainsContainer:
    def __init__(self):
        self.map = {}
        self.uidmap = []
        self.tmpuid = ""
        self.position = -1

    def add_values(self, ipr, uid, start):
        if self.tmpuid != uid:
            self.position += 1
            self.uidmap.insert(self.position, uid)
            self.tmpuid = uid

        if ipr in self.map:
            ht = self.map[ipr]
            if self.position in ht:
                ht[self.position].append(start)
            else:
                starts = [start]
                ht[self.position] = starts
        else:
            ht = {}
            starts = [start]
            ht[self.position] = starts
            self.map[ipr] = ht

    def get_uid_with_ipr(self, ipr):
        ht = {}
        iprdata = self.map.get(ipr)

        for i in iprdata.keys():
            ht[self.uidmap[i]] = iprdata[i]

        return ht

    def get_uid(self, i):
        return self.uidmap[i]

    def save(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)


class Domains:
    def __init__(self):
        print('initialising domains')
        self.args = OptionsParser()
        self.setting = Settings().init()

        if not OptionsParser().parse_options(self.setting):
            sys.exit(1)

        self.setting.set_work_dir()
        create_directory(self.setting.get_work_dir())

        mongo_client = MongoClient(self.setting.get_host())
        self.database = mongo_client.get_database(self.setting.get_db())

    def run_domains(self):
        print('running domains')
        cores = multiprocessing.cpu_count()
        if self.setting.get_threads() > cores:
            raise ValueError(f"Exceeded the maximum number of available processors. Available processors: {cores}.")
        self.setting.set_num_threads_for_external_process(cores)
        os.environ["JAVA_OPTS"] = f"-Djava.util.concurrent.ForkJoinPool.common.parallelism={self.setting.get_threads()}"
        templates_for_thread_list = []
        blast_prot_map = {}
        thread_list = []

        domains_map = read_dom_file(self.setting.get_dom_file())
        templates = list(domains_map.keys())

        if len(self.setting.get_blast_file()) > 0:
            blast_prot_map = parse_blast_file(self.setting.get_blast_file(), self.setting)

        for i, template in enumerate(templates):
            list_num = i % self.setting.get_threads()

            if list_num == len(templates_for_thread_list):
                templates_for_thread_list.append([])

            if templates_for_thread_list[list_num] is None:
                templates_for_thread_list[list_num] = []

            if len(template) >= 5:
                templates_for_thread_list[list_num].append(template)

        for i in range(len(templates_for_thread_list)):
            thread_list.append(DomainsThread(domains_map, blast_prot_map, templates_for_thread_list[i], self.database,
                                             self.setting.get_nwth()))

        results = merge_thread_results(thread_list)

        self.write_results(results, self.database.get_collection(self.setting.get_collection_goa()))

        self.setting.get_out().close()

    def write_results(self, results, collgoa):
        for template in results.keys():
            self.setting.get_out().write(f">{template}\n")

            for uid in results[template].get_neighbor_map().keys():
                if results[template].get_prot_custers() is not None and uid in results[template].get_prot_custers():
                    for go_id in results[template].get_prot_custers()[uid].get_gos():
                        self.setting.get_out().write(
                            f"{go_id}\t{results[template].get_neighbor_map()[uid]}\t{uid}\t0\n")
                else:
                    for doc in collgoa.find({"uid": uid}):
                        go_id = doc["goid"]

                        if uid not in results[template].get_prot_custers() or go_id not in \
                                results[template].get_prot_custers()[uid].get_gos():
                            self.setting.get_out().write(
                                go_id + "\t" + str(results[template].get_neighbor_map()[uid]) + "\t" + uid + "\t0\n")
            self.setting.get_out().flush()

    def get_settings(self):
        return self.setting


class DomainsThread(Process):
    def __init__(self, domains_map, blast_prot_map, templates, database, nwth):
        super().__init__()
        self.domains_map = domains_map
        self.blast_prot_map = blast_prot_map
        self.templates = templates
        self.database = database
        self.nwth = nwth
        self.output = Manager().dict()

    def run(self):
        for template in self.templates:
            if self.blast_prot_map.get(template) is not None:
                self.get_neighbors(self.domains_map.get(template), template, self.blast_prot_map.get(template))

    def get_neighbors(self, mapp, pid, blast_prot):
        try:
            doms = self.get_domains(mapp)
            dd = DataDomains(self.database, doms)
            if not dd.get_exit():
                uidscore = dd.get_results().get_uid_over(self.nwth)
                if uidscore:
                    at = ArgotGraph(uidscore, pid, self.database, blast_prot, self.ident)
                    if pid not in self.output:
                        self.output[pid] = ThreadData(at.get_prot_clusters(), at.get_cluster())
            else:
                logging.error(f"ERROR: suck: {pid}")
        except Exception as e:
            logging.error(f"ERROR: {e}, {pid}")

    @staticmethod
    def get_domains(mapp):
        buffer = []
        for p in sorted(mapp.keys()):
            if len(mapp[p]) > 1:
                for ipr in sorted(list(mapp[p])):
                    buffer.append(ipr)
            else:
                buffer.extend(list(mapp[p]))
        return buffer


class Edge:
    def __init__(self, node, rel):
        self.node = node
        self.rel = rel

    def get_relationship_type(self):
        return self.rel

    def __str__(self):
        return str(self.rel)

    def get_node(self):
        return self.node


class GOGraphOWL:
    def __init__(self, owl_file, g, e):
        self.typeg = g
        self.typee = e
        self.node_map = {}
        self.relation = Relation().instance()
        self.load_from_file(owl_file)
        self.model = None
        self.descendantb = None
        self.ancestorb = None
        self.visited = set()
        self.id = 0
        self.gonode = None

    def load_from_file(self, file):
        self.model = Graph()
        self.model.parse(file, format="xml")
        for subject, predicate, obj in self.model.triples((None, OWL.Class, None)):
            self.browse_class(subject)

    def browse_class(self, ontclass):
        localname = ontclass.get_local_name()
        localname = localname.replace('_', ':')
        if localname.startswith("GO:"):
            ontpr = self.model.getOntProperty("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace")
            namespace = ontclass.get_property_value(ontpr)
            if namespace is not None:
                if localname in self.node_map:
                    self.gonode = self.node_map[localname]
                else:
                    self.gonode = self.create_gonode(localname)
                self.gonode.setNamespace(str(namespace))
                ontpr = self.model.getOntProperty("http://www.w3.org/2000/01/rdf-schema#label")
                label = ontclass.get_property_value(ontpr)
                self.gonode.set_comment(str(label))
                iterst = ontclass.listProperties()

                while iterst.hasNext():
                    stmt = iterst.nextStatement()
                    predicate = stmt.getPredicate()
                    if predicate.hasURI("http://www.w3.org/2000/01/rdf-schema#subClassOf"):
                        rdfnode = stmt.getObject()
                        if rdfnode.isURIResource():
                            ontid = stmt.getObject().asResource().getLocalName()
                            ontid = ontid.replace('_', ':')
                            self.add_node(ontid, RelationshipType.is_a)
                        else:
                            oc = rdfnode.asOntClass()
                            if oc.isRestriction():
                                r = oc.asRestriction()
                                if r.isSomeValuesFromRestriction():
                                    av = r.asSomeValuesFromRestriction()
                                    ontid = av.getSomeValuesFrom().getLocalName()
                                    ontid = ontid.replace('_', ':')
                                    if av.getOnProperty().getLocalName() in Relation:
                                        rel = RelationshipType[av.getOnProperty().getLocalName()]
                                        self.add_node(ontid, rel)
                    elif predicate.hasURI("http://www.w3.org/2002/07/owl#equivalentClass"):
                        itr = ontclass.listEquivalentClasses()
                        while itr.hasNext():
                            ontc = itr.next()
                            inc = ontc.asIntersectionClass()
                            rdflist = inc.getOperands()
                            for ii in rdflist.iterator():
                                n = ii.next()
                                if n.isResource():
                                    res = n.asResource()
                                    oc = res.asOntClass()
                                    if oc.isRestriction():
                                        r = oc.asRestriction()
                                        if r.isSomeValuesFromRestriction():
                                            av = r.asSomeValuesFromRestriction()
                                            ontid = av.getSomeValuesFrom().getLocalName()
                                            ontid = ontid.replace('_', ':')
                                            if av.getOnProperty().getLocalName() in Relation:
                                                rel = RelationshipType[av.getOnProperty().getLocalName()]
                                                self.add_node(ontid, rel)
                    elif predicate.hasURI("http://www.w3.org/2002/07/owl#disjointWith"):
                        ontid = stmt.getObject().asResource().getLocalName()
                        if ontid is not None:
                            ontid = ontid.replace('_', ':')
                            if ontid in self.node_map:
                                n = self.node_map[ontid]
                                n.addDisjoint(self.gonode.get_ont_id())
                            else:
                                assert ontid.startswith("GO:")
                                self.create_gonode(ontid.addDisjoint(self.gonode.get_ont_id()))
                            self.gonode.addDisjoint(ontid)
        elif localname.startswith("CHEBI:") or localname.startswith("PR:"):
            if localname in self.node_map:
                extnode = self.node_map[localname]
            else:
                extnode = NodeFactory().get_node(self.typee)
                extnode.set_id(self.id + 1)
                ns = localname.split(":")
                extnode.set_name_space(ns[0])
                self.node_map[localname] = extnode
            extnode.set_ont_id(localname)

    def add_node(self, ontid, rel):
        if ontid.startswith("GO:"):
            if ontid in self.node_map:
                gn = self.node_map[ontid]
            else:
                gn = self.create_gonode(ontid)
            self.set_node(gn, self.gonode, rel)
        elif ontid.startswith("CHEBI:") or ontid.startswith("PR:"):
            if ontid in self.node_map:
                extn = self.node_map[ontid]
            else:
                extn = NodeFactory().get_node(self.typee)
                extn.set_ont_id(ontid)
                self.id += 1
                extn.set_id(self.id)
                ns = ontid.split(":")
                extn.set_namespace(ns[0])
                self.node_map[ontid] = extn
            self.set_node(self.gonode, extn, rel)

    @staticmethod
    def set_node(parent, child, rel):
        if parent.contains_child(child) and child.contains_parent(parent):
            return

        edgep = Edge(parent, rel)
        edgec = Edge(child, rel)
        parent.add_child(edgec)
        child.add_parent(edgep)
        parent.add_child_node(child, edgec)
        child.add_parent_node(parent, edgep)

    def create_gonode(self, goid):
        n = NodeFactory().get_node(self.typeg)
        n.set_ont_id(goid)
        self.id += 1
        n.set_id(self.id)
        self.node_map[goid] = n
        return n

    def get_node_map(self):
        return self.node_map

    @staticmethod
    def get_parents(node):
        return node.get_parents()

    @staticmethod
    def get_children(node):
        return node.get_children()

    @staticmethod
    def get_go_children(node):
        listt = []
        for e in node.get_children():
            if Relation().get_relation(e.get_relationship_type()) != Direction.propagationUp:
                if isinstance(e.get_node(), Node):
                    listt.append(e)
        return listt

    @staticmethod
    def get_go_parents(node):
        listt = []
        for e in node.get_parents():
            if Relation().get_relation(e.get_relationship_type()) != Direction.propagationDown:
                if isinstance(e.get_node(), Node):
                    listt.append(e)
        return listt

    def is_ancestor_of(self, ancestor, node):
        self.browse_ancestor(ancestor, node)
        if self.ancestorb:
            self.ancestorb = False
            return True
        return False

    def browse_ancestor(self, ancestor, node):
        for e in node.get_parents():
            if Relation().get_relation(e.get_relationship_type()) != Direction.propagationDown:
                if e.get_node() == ancestor:
                    self.ancestorb = True
                    return
                self.browse_ancestor(ancestor, e.get_node())

    def is_descendant_of(self, descendant, node):
        self.browse_descendant(descendant, node)
        self.visited.clear()
        if self.descendantb:
            self.descendantb = False
            return True
        return False

    def browse_descendant(self, descendant, node):
        if node not in self.visited:
            self.visited.add(node)
            for e in node.get_children():
                if Relation().get_relation(e.get_relationship_type()) != Direction.propagationUp:
                    if e.get_node() == descendant:
                        self.descendantb = True
                        return
                    self.browse_descendant(descendant, e.get_node())

    @staticmethod
    def process_id():
        return "GO:0008150"

    @staticmethod
    def function_id():
        return "GO:0003674"

    @staticmethod
    def component_id():
        return "GO:0005575"

    def clean_all(self):
        for key in self.node_map.keys():
            n = self.node_map[key]
            if n.is_flagged():
                n.clean()

    def read_disjoint(self, disj_file):
        try:
            with open(disj_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    line = line.strip()
                    if self.exists(line):
                        try:
                            n = self.get_go_node(line)
                            n.set_int_disjoint(True)
                        except NodeNotFoundException as ex:
                            logging.getLogger(GOGraphOWL.__name__).log(logging.CRITICAL, None, ex)
        except (FileNotFoundError, IOError) as ex:
            logging.getLogger(GOGraphOWL.__name__).log(logging.CRITICAL, None, ex)

    def dump_sif(self):
        buffer = set()
        try:
            with open("map.sif", "w") as fw:
                for ontid in self.node_map.keys():
                    n = self.node_map[ontid]
                    buffer.add(ontid)
                    for edge in n.get_parents():
                        if edge.get_node().get_ont_id() not in buffer:
                            fw.write(f"{ontid}\t{edge.get_relationship_type()}\t{edge.get_node().get_ont_id()}\n")
                    for edge in n.get_children():
                        if edge.get_node().get_ont_id() not in buffer:
                            fw.write(f"{edge.get_node().get_ont_id()}\t{edge.get_relationship_type()}\t{ontid}\n")
        except IOError as ex:
            logging.getLogger(GOGraphOWL.__name__).log(logging.CRITICAL, None, ex)

    def dump_edge(self):
        try:
            with open("dataEdge.txt", "w") as fw:
                for ontid in self.node_map.keys():
                    n = self.node_map[ontid]
                    if isinstance(n, Node) and len(n.get_parents()) > 1:
                        fw.write(n.get_ont_id())
                        for edge in n.get_parents():
                            fw.write(f"\t{edge.get_node().get_ont_id()}")
                        fw.write("\n")
        except IOError as ex:
            logging.getLogger(GOGraphOWL.__name__).log(logging.CRITICAL, None, ex)

    def remove_cross_path(self):
        for ontid in self.node_map.keys():
            if self.node_map[ontid].is_go_node():
                node = self.node_map[ontid]
                buffer = []
                for e in node.get_parents():
                    if e.get_node().is_go_node():
                        n = e.get_node()
                        if node.get_name_space() != n.get_name_space():
                            buffer.append(e)
                            n.remove_child_node(node)
                            print(f"{node.get_ont_id()} -- {n.get_ont_id()}")
                for e in buffer:
                    node.remove_parent(e)

    def get_node(self, ontid):
        if ontid in self.node_map:
            n = self.node_map[ontid]
        else:
            raise NodeNotFoundException("The Ontology ID " + ontid + " is not present in the current ontology.")
        return n

    def get_go_node(self, goid):
        if goid in self.node_map:
            n = self.node_map[goid]
            if isinstance(n, Node):
                return n
        else:
            raise NodeNotFoundException("The GO ID " + goid + " is not present in the current ontology.")
        return None

    def get_ancestors_of(self, node):
        anc = []
        self.browse_ancestors(anc, node)
        return anc

    def browse_ancestors(self, anc, node):
        if not node.is_root():
            for e in node.get_parents():
                if self.relation.getRelation(e.get_relationship_type()) != Direction.propagationDown:
                    if isinstance(e.get_node(), Node):
                        anc.append(e.get_node())
                    self.browse_ancestors(anc, e.get_node())

    def get_all_go_descendants(self, n):
        listt = set()
        if not n.is_leaf():
            for e in n.get_children():
                if isinstance(e.get_node(), Node):
                    self.walk_descendants(e.get_node(), listt)
        return listt

    def walk_descendants(self, n, listt):
        listt.add(n)
        if not n.is_leaf():
            for e in n.get_children():
                if isinstance(e.get_node(), Node):
                    if self.relation.get_relation(e.get_relationship_type()) != Direction.propagationUp:
                        self.walk_descendants(e.get_node(), listt)

    def get_all_go_ancestors(self, n):
        listt = set()
        if not n.is_root():
            for e in n.get_parents():
                if isinstance(e.get_node(), Node):
                    self.walk_ancestors(e.get_node(), listt)
        return listt

    def walk_ancestors(self, n, listt):
        listt.add(n)
        if not n.is_root():
            for e in n.get_parents():
                if isinstance(e.get_node(), Node):
                    if self.relation.get_relation(e.get_relationship_type()) != Direction.propagationDown:
                        self.walk_ancestors(e.get_node(), listt)

    def exists(self, ontid):
        return ontid in self.node_map

    def get_root(self, namespace):
        gid = ""
        root = None
        if namespace == OntNamespace.molecular_function:
            gid = self.function_id()
        elif namespace == OntNamespace.biological_process:
            gid = self.process_id()
        elif namespace == OntNamespace.cellular_component:
            gid = self.component_id()

        try:
            root = self.node_map[gid]
        except NodeNotFoundException as ex:
            logging.getLogger(GOGraphOWL.__name__).log(logging.CRITICAL, None, ex)
        return root


class GraphFreq:
    _instance = None

    def __init__(self):
        self.total_terms = None
        self.id2freq = {}
        self.graph = None
        self.max_bp_ic = 0.0
        self.max_mf_ic = 0.0
        self.max_cc_ic = 0.0
        self.relation = Relation()
        self._instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls.relation = Relation.instance()
            cls.id2freq = {}
        return cls._instance

    @classmethod
    def instance(cls):
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def init(self, graph, database, collection):
        self.graph = graph
        try:
            self.load_freq_from_mongo(database, collection)
            for goid, freq in self.id2freq.items():
                if graph.exists(goid):
                    self.add_abs_value(graph.getGONode(goid), freq)

            self.calc_freq()
        except NodeNotFoundException as ex:
            print('o forse era qui?')

    def add_abs_value(self, node, w):
        node.addAbsValue(w)
        for e in self.graph.get_go_parents(node):
            if self.relation.get_relation(e.getRelationshipType()) != Direction.propagationDown:
                self.add_abs_value(e.get_node(), w)

    def calc_freq(self):
        mapnode = self.graph.get_node_map()
        maxw = max(self.graph.get_root(OntNamespace.biological_process).get_abs_value(),
                   self.graph.get_root(OntNamespace.molecular_function).get_abs_value())
        maxw = max(maxw, self.graph.get_root(OntNamespace.cellular_component).get_abs_value())
        for ontid, node in mapnode.items():
            if isinstance(node, Node):
                if node.get_ont_id() in self.id2freq:
                    node.set_real_count(self.id2freq[node.get_ont_id()])
                    absfreq = node.get_real_count() / self.total_terms
                    node.set_real_frequency(absfreq)

                freq = node.get_abs_value() / maxw

                if freq > 0:
                    ic = -log(freq)
                    if node.get_namespace() == OntNamespace.biological_process:
                        if ic > self.max_bp_ic:
                            self.max_bp_ic = ic
                    elif node.get_namespace() == OntNamespace.molecular_function:
                        if ic > self.max_mf_ic:
                            self.max_mf_ic = ic
                    elif node.get_namespace() == OntNamespace.cellular_component:
                        if ic > self.max_cc_ic:
                            self.max_cc_ic = ic
                node.set_frequency(freq)

    def load_freq_from_mongo(self, database, collection):
        coll = database[collection]
        cursor = coll.find(projection={'_id': False})
        for doc in cursor:
            t = doc["freq"]
            self.id2freq[doc["goid"]] = t
            self.total_terms += t

    def max_bp_ic(self):
        return self.max_bp_ic

    def max_mf_ic(self):
        return self.max_mf_ic

    def max_cc_ic(self):
        return self.max_cc_ic

    def max_ic(self, namespace):
        if namespace == OntNamespace.biological_process:
            return self.max_bp_ic
        elif namespace == OntNamespace.molecular_function:
            return self.max_mf_ic
        elif namespace == OntNamespace.cellular_component:
            return self.max_cc_ic
        else:
            return 0

    def get_total_term_number(self):
        return self.total_terms


class LoadInRam:
    def __init__(self, filename=None):
        self.baos = io.BytesIO()
        if filename:
            self.write_bytes(filename)

    def write_bytes(self, data):
        if isinstance(data, bytes):
            self.baos.write(data)
        elif isinstance(data, str):
            with open(data, 'rb') as f:
                self.baos.write(f.read())

    def get_content(self):
        return self.baos.getvalue().decode()

    def get_byte_array(self):
        return self.baos.getvalue()

    def size(self):
        return self.baos.getbuffer().nbytes

    def reset(self):
        self.baos.truncate(0)
        self.baos.seek(0)

    def close(self):
        self.baos.close()

    def write_to_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.get_content())


class Igraph:
    def fast_greedy(self, edges, res, lenn):
        pass


if os.name == 'nt':  # Windows
    libgraph = CDLL('GraphMod.dll')
elif os.name == 'posix':  # Linux, macOS
    libgraph = CDLL('libGraphMod.so')
else:
    raise OSError('Unsupported operating system')
libgraph.fastGreedy.argtypes = [POINTER(c_int), POINTER(c_int), c_int]
libgraph.fastGreedy.restype = None


class MakeClusters:
    def get_fast_greedy_res(self, edges, reslen):
        # Set the library path and load the Igraph interface
        os.environ['PATH'] = os.environ.get('PATH', '') + os.pathsep + os.path.abspath('libgraph')
        igraph = Igraph()

        # Call the fastGreedy function and return the result
        res = (c_int * reslen)()
        igraph.fast_greedy(edges, res, len(edges))
        return res


class Modularity:
    def __init__(self, graph, c):
        self.map = {}
        self.buffer = []
        edges = np.zeros(graph.edge_set().size() * 2, dtype=int)
        index = np.zeros(graph.vertex_set().size(), dtype=c)

        i = 0
        n = 0
        edgeset = graph.edge_set()
        for e in edgeset:
            node = graph.get_edge_source(e)
            if node not in self.map:
                index[i] = node
                self.map[node] = i
                i += 1
            edges[n] = self.map[node]
            n += 1
            node = graph.get_edge_target(e)
            if node not in self.map:
                index[i] = node
                self.map[node] = i
                i += 1
            edges[n] = self.map[node]
            n += 1
        self.map = {}  # deleting
        ss = MakeClusters()
        cluster = ss.get_fast_greedy_res(edges, index.size)
        memberships = {}
        for i in range(cluster.__len__()):
            if cluster[i] in memberships:
                memberships[cluster[i]].add(index[i])
            else:
                hs = set()
                hs.add(index[i])
                memberships[cluster[i]] = hs
        alledges = graph.edge_set()

        for hs in memberships.values():
            g = SimpleWeightedGraph(DefaultWeightedEdge)
            for e in alledges:
                k = graph.get_edge_source(e)
                m = graph.get_edge_target(e)
                if k in hs and m in hs:
                    if k not in g:
                        g.add_vertex(k)
                    if m not in g:
                        g.add_vertex(m)
                    dwe = g.add_edge(k, m)
                    g.set_edge_weight(dwe, graph.get_edge_weight(e))
            self.buffer.append(g)

    def get_clusters(self):
        return self.buffer


class NodeFactory:
    @staticmethod
    def get_node(node):
        try:
            instance = node()
            return instance
        except Exception as ex:
            logging.getLogger(NodeFactory.__name__).log(logging.ERROR, ex)
        return None


class NodeNotFoundException(Exception):
    def __init__(self, msg=None):
        super().__init__(msg)


class NW:
    def __init__(self, x, y, weightscore):
        setting = Settings().init()
        self.x = x
        self.y = y
        self.weight_score = weightscore
        self.gap_score = setting.get_gap_score()
        self.match_score = setting.get_match_score()
        self.mismatch_score = setting.get_mismatch_score()
        self.x_len = len(x)
        self.y_len = len(y)
        self.score_array = [[0.0] * (self.x_len + 1) for _ in range(self.y_len + 1)]
        self.fill_score_array()

    def fill_score_array(self):
        # Fill the top row and left column:
        for col in range(self.x_len + 1):
            self.score_array[0][col] = self.gap_score * col
        for row in range(self.y_len + 1):
            self.score_array[row][0] = self.gap_score * row

        # Now fill in the rest of the array:
        for row in range(1, self.y_len + 1):
            for col in range(1, self.x_len + 1):
                if self.x[col - 1] == self.y[row - 1]:
                    northwest = self.score_array[row - 1][col - 1] + self.match_score * (
                                2 * self.weight_score[self.x[col - 1]])
                else:
                    northwest = self.score_array[row - 1][col - 1] + self.mismatch_score * (
                                self.weight_score[self.y[row - 1]] + self.weight_score[self.x[col - 1]])
                west = self.score_array[row][col - 1] + self.gap_score
                north = self.score_array[row - 1][col] + self.gap_score
                best = northwest
                if north > best:
                    best = north
                if west > best:
                    best = west
                self.score_array[row][col] = best

    def get_similarity(self):
        sim = 0.0
        if self.y_len == 0 and self.x_len == 0:
            return sim

        score = self.score_array[self.y_len][self.x_len]
        smn = self.x_len * self.gap_score + self.y_len * self.gap_score
        selfx = 0.0
        for dom in self.x:
            selfx += self.match_score * self.weight_score[dom] * 2
        selfy = 0.0
        for dom in self.y:
            selfy += self.match_score * self.weight_score[dom] * 2
        sim = min((score - smn) / (selfx - smn), (score - smn) / (selfy - smn))
        return sim


class OntNamespace(Enum):
    molecular_function = 1
    biological_process = 2
    cellular_component = 3
    CHEBI = 4
    PR = 5
    GOCHE = 6
    BFO = 7
    OBA = 8
    CL = 9
    SO = 10
    UBERON = 11
    PO = 12
    PATO = 13
    NCBITaxon = 14


class OptionsParser:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.parser = argparse.ArgumentParser(description='Description of your program')
        self.parser.add_argument('-i', '--input', help='The output of interproscan (mandatory)')
        self.parser.add_argument('-s', '--setting', help='The settings file (mandatory)')
        self.parser.add_argument('-b', '--blast', help='The output of Blast', nargs='?')
        self.parser.add_argument('-w', '--nwth', help='The score alignment to use as threshold')
        self.parser.add_argument('-o', '--output', help='The output file')
        self.parser.add_argument('-d', '--seed_dom',
                                 help='The algorithm to use to individuate the seed domains <minfreq (default), maxfreq, all>')
        self.parser.add_argument('-p', '--print', help='Print gml graph', action='store_true')
        self.parser.add_argument('-r', '--rom', help='Use ROM for temporary files instead of RAM', nargs='?')
        self.parser.add_argument('-t', '--threads', help='Number of threads', nargs='?')
        self.args = self.parser.parse_args()

    def parse_options(self, settings):
        if self.args.input:
            settings.set_dom_file(self.args.input)
        else:
            self.logger.error('Input file is missing')
            return False

        if self.args.setting:
            settings.read_settings(self.args.setting)
        else:
            if os.path.exists(settings.FILE):
                settings.read_settings()
            else:
                self.logger.error('Setting file is missing')
                return False

        if self.args.nwth:
            settings.set_nwth(float(self.args.nwth))

        if self.args.output:
            try:
                settings.set_out(open(self.args.output, 'w'))
            except IOError as e:
                self.logger.error('Output file cannot be opened: {}'.format(e))
                self.parser.print_help()
                return False
        else:
            settings.setOut(sys.stdout)

        if self.args.seed_dom:
            alg = self.args.seed_dom
            if alg in ['minfreq', 'maxfreq', 'all']:
                settings.set_seed_dom(alg)
            else:
                self.logger.error('Invalid seed domain algorithm specified')
                self.parser.print_help()
                return False
        else:
            settings.set_seed_dom('minfreq')

        if self.args.print:
            settings.set_print_gml(True)

        if self.args.rom:
            settings.use_rom()

        if self.args.blast:
            settings.set_blast_file(self.args.blast)

        if self.args.threads:
            settings.set_threads(int(self.args.threads))

        return True


class Direction(Enum):
    propagationUp = 1
    propagationDown = 2
    transitive = 3
    symmetric = 4


class RelationshipType(Enum):
    is_a = 1
    BFO_0000050 = 2
    BFO_0000051 = 3
    BFO_0000066 = 4
    RO_0002211 = 5
    RO_0002212 = 6
    RO_0002213 = 7
    RO_0002215 = 8
    RO_0002216 = 9
    RO_0002091 = 10
    RO_0002092 = 11
    RO_0002093 = 12


class Relation:
    relation = None
    rel = None
    relations = None

    def __init__(self):
        self.relation = {}
        self.relations = set()
        self.relation[RelationshipType.is_a] = Direction.transitive
        self.relation[RelationshipType.BFO_0000050] = Direction.transitive
        self.relation[RelationshipType.BFO_0000051] = Direction.propagationDown  # non usare has_part
        self.relation[RelationshipType.BFO_0000066] = Direction.transitive
        self.relation[RelationshipType.RO_0002211] = Direction.transitive
        self.relation[RelationshipType.RO_0002212] = Direction.transitive
        self.relation[RelationshipType.RO_0002213] = Direction.transitive
        self.relation[RelationshipType.RO_0002215] = Direction.transitive
        self.relation[RelationshipType.RO_0002216] = Direction.transitive
        self.relation[RelationshipType.RO_0002091] = Direction.transitive
        self.relation[RelationshipType.RO_0002092] = Direction.transitive
        self.relation[RelationshipType.RO_0002093] = Direction.transitive
        self.relations.add("is_a")
        self.relations.add("BFO_0000050")  # part_of
        self.relations.add("BFO_0000066")  # occurs_in
        self.relations.add("RO_0002211")  # regulates
        self.relations.add("RO_0002212")  # negatively regulates
        self.relations.add("RO_0002213")  # positively regulates
        self.relations.add("RO_0002215")  # capable of
        self.relations.add("RO_0002216")  # capable of part of
        self.relations.add("RO_0002091")
        self.relations.add("RO_0002092")
        self.relations.add("RO_0002093")

    @staticmethod
    def instance():
        if not Relation.rel:
            Relation.rel = Relation()
        return Relation.rel

    def get_relation(self, reltype):
        return self.relation[reltype]

    def contains(self, rel):
        return rel in self.relations


class RetrieveProteins:
    instance = None

    def __init__(self, mongo=False):
        self.inputmap = {}
        self.offset = {}
        self.settings = Settings()
        self.db_input = self.settings.get_input_fasta()
        self.db_file = self.settings.get_fasta_db()
        self.read_index()
        if not mongo:
            self.read_uniprot_index()

    @staticmethod
    def init(mongo=False):
        if RetrieveProteins.instance is None:
            RetrieveProteins.instance = RetrieveProteins(mongo)
        return RetrieveProteins.instance

    def read_index(self):
        index_file_input = self.db_input + ".fai"
        with open(index_file_input, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                data = line.split("\t")
                pid = data[0]
                self.inputmap[pid] = int(data[1])

    def read_uniprot_index(self):
        index_file = self.db_file + ".fai"
        with open(index_file, "r") as fr:
            for line in fr:
                if line.startswith("#"):
                    continue
                line = line.strip()
                data = line.split("\t")
                self.offset[data[0]] = int(data[1])

    def write_query(self, rafinput, output, pid):
        lenseq = 0
        if pid not in self.inputmap:
            print("ERROR: no such protein in index: " + pid, file=sys.stderr)
            sys.exit(1)
        offs = self.inputmap[pid]
        rafinput.seek(offs)
        line = rafinput.readline().decode().strip()
        output.write(">" + pid + "\n")
        while line:
            line = line.strip()
            if line.startswith(">"):
                break
            elif len(line) == 0:
                continue
            lenseq += len(line)
            output.write(line + "\n")
            line = rafinput.readline().decode().strip()

        rafinput.seek(self.inputmap[pid])
        line = rafinput.readline()
        output.write(line.decode('utf-8'))

        return lenseq + len(line)

    def write_target(self, raf, output, uid, lenn, offs):
        minlen = 0
        maxlen = 0
        perc = self.settings.get_target_perc_len()
        sb = ""
        if perc > 0:
            minlen = int(lenn - lenn * perc)
            maxlen = int(lenn + lenn * perc)

        if offs != -1:
            raf.seek(offs)
            line = raf.readline()
            if line:
                while not line.startswith(">"):
                    line = raf.readline()
                    if not line:
                        return False
                sb += line.decode('utf-8').strip() + "\n"
                len_line = len(line)
                while len(line) == len_line:
                    line = raf.readline()
                    sb += line.decode('utf-8').strip() + "\n"
                if minlen == 0 or minlen <= len(sb) <= maxlen:
                    output.write(">" + uid + "\n")
                    output.write(sb)
                    return True

        return False

    def write_prot_file(self, uids, outfile, pid, database):
        count = 1
        with open(outfile, 'w') as output, open(self.db_file, 'rb') as raf, open(self.db_input, 'rb') as rafinput:
            len_line = self.write_query(rafinput, output, pid)
            collection = database[Settings().init().get_collection_uniprotx()]
            for uid in uids:
                cursor = collection.find({"uid": uid})
                if cursor.count() > 0:
                    doc = cursor.next()
                    if doc:
                        offs = -1
                        res = doc.get("index")
                        if isinstance(res, int):
                            offs = res
                        if count <= self.settings.get_maxseqs():
                            if self.write_target(raf, output, uid, len_line, offs):
                                count += 1
                        else:
                            break

    def write_prot_file_repr(self, uids, outfile, pid, database):
        count = 1

        with open(self.db_file, 'rb') as raf:
            with open(self.db_input, 'rb') as rafinput:
                output = open(outfile, 'w')
                lenn = self.write_query(rafinput, output, pid)
                collection = database.get_collection(Settings().init().get_collection_uniprotx())
                for uid in uids:
                    cursor = collection.find({"uid": uid}).limit(1)
                    if cursor.count() > 0:
                        doc = cursor[0]
                        if not doc.empty:
                            offs = -1
                            res = doc.get("index")
                            if isinstance(res, int):
                                offs = res
                            elif isinstance(res, float):
                                offs = int(res)

                            if count <= self.settings.get_max_num_repr():
                                if self.write_target(raf, output, uid, lenn, offs):
                                    count += 1
                            else:
                                break

                output.close()

    def write_single_prot(self, uid, outfile, database):
        with open(self.db_file, 'rb') as raf:
            output = open(outfile, 'w')
            collection = database.get_collection(Settings().init().get_collection_uniprotx())
            cursor = collection.find({"uid": uid}).limit(1)
            if cursor.count() > 0:
                doc = cursor[0]
                if not doc.empty:
                    offs = -1
                    res = doc.get("index")
                    if isinstance(res, int):
                        offs = res
                    elif isinstance(res, float):
                        offs = int(res)

                    self.write_target(raf, output, uid, 0, offs)

            output.close()


class Settings:
    settings = None
    FILE = "setting.conf"

    def __init__(self):
        self.gl_search_align_cmd = None
        self.mm_search_align_cmd = None
        self.fasta_db = None
        self.work_dir = None
        self.input_fasta = None
        self.host = None
        self.db = None
        self.coll_freq = None
        self.coll_interpro = None
        self.coll_goa = None
        self.uniprotx = None
        self.seed_dom = None
        self.max_iter = 2
        self.match_score = None
        self.gap_score = None
        self.mismatch_score = None
        self.nwth = None
        self.perc_len = None
        self.sim_thr = None
        self.max_seqs = 0
        self.deep_thr = 50
        self.print_gml = False
        self.use_rom = False
        self.hit_num = None
        self.min_score = None
        self.max_mum_repr = 1000
        self.minum_go = None
        self.blast_thr = 90
        self.threads = 1
        self.dom_file = None
        self.blast_file = ""
        self.out = None
        self.num_threads_for_external_process = 1

    @staticmethod
    def init():
        if Settings.settings is None:
            Settings.settings = Settings()
        return Settings.settings

    def read_settings(self, file=FILE):
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                cmd, value = line.split('=', 1)
                cmd = cmd.strip()
                value = value.strip()
                if cmd == "GLSEARCHALIGN":
                    self.gl_search_align_cmd = value
                elif cmd == "MMSEARCHALIGN":
                    self.mm_search_align_cmd = value
                elif cmd == "FASTADB":
                    self.fasta_db = value
                elif cmd == "WORKDIR":
                    self.work_dir = value
                elif cmd == "INPUTFASTA":
                    self.input_fasta = value
                elif cmd == "HOST":
                    self.host = value
                elif cmd == "DB":
                    self.db = value
                elif cmd == "PERCLEN":
                    self.perc_len = float(value)
                elif cmd == "SIMTHR":
                    self.sim_thr = float(value)
                elif cmd == "MAXITER":
                    self.max_iter = int(value)
                elif cmd == "MATCHSCORE":
                    self.match_score = float(value)
                elif cmd == "GAPSCORE":
                    self.gap_score = float(value)
                elif cmd == "MISMATCHSCORE":
                    self.mismatch_score = float(value)
                elif cmd == "NWTH":
                    self.nwth = float(value)
                elif cmd == "FREQ":
                    self.coll_freq = value
                elif cmd == "INTERPRO":
                    self.coll_interpro = value
                elif cmd == "GOA":
                    self.coll_goa = value
                elif cmd == "UNIPROTX":
                    self.uniprotx = value
                elif cmd == "MAXSEQS":
                    self.max_seqs = int(value)
                elif cmd == "DEEPTHR":
                    self.deep_thr = float(value)
                elif cmd == "HITNUM":
                    self.hit_num = int(value)
                elif cmd == "MINSCORE":
                    self.min_score = float(value)
                elif cmd == "MAXNUMREPR":
                    self.max_num_repr = int(value)
                    if self.max_seqs == 0:
                        self.max_seqs = 10 * self.max_num_repr
                elif cmd == "MINUMGO":
                    self.minum_go = int(value)
                elif cmd == "BLASTTHR":
                    self.blast_thr = float(value)
        if self.max_seqs == 0:
            self.max_seqs = 4 * self.max_num_repr

    def get_glsearch_alignment_cmd(self):
        return self.gl_search_align_cmd

    def get_mmsearch_alignment_cmd(self):
        return self.mm_search_align_cmd

    def get_fasta_db(self):
        return self.fasta_db

    def get_work_dir(self):
        return self.work_dir

    def get_input_fasta(self):
        return self.input_fasta

    def get_host(self):
        return self.host

    def get_db(self):
        return self.db

    def get_nwth(self):
        return self.nwth

    def set_nwth(self, nwth):
        self.nwth = nwth

    def get_collection_frequences(self):
        return self.coll_freq

    def get_collection_interpro(self):
        return self.coll_interpro

    def get_collection_goa(self):
        return self.coll_goa

    def get_collection_uniprotx(self):
        return self.uniprotx

    def get_seed_dom(self):
        return self.seed_dom

    def set_seed_dom(self, seed_dom):
        self.seed_dom = seed_dom

    def get_match_score(self):
        return self.match_score

    def get_gap_score(self):
        return self.gap_score

    def get_mismatch_score(self):
        return self.mismatch_score

    def get_target_perc_len(self):
        return self.perc_len

    def get_max_iteration(self):
        return self.max_iter

    def get_similarity_thr(self):
        return self.sim_thr

    def get_maxseqs(self):
        return self.max_seqs

    def get_deepthr(self):
        return self.deep_thr

    def get_hitnum(self):
        return self.hit_num

    def get_minscore(self):
        return self.min_score

    def get_minumgo(self):
        return self.minum_go

    def is_rom_used(self):
        return self.use_rom

    def is_printgml(self):
        return self.print_gml

    def set_printgml(self, printgml):
        self.print_gml = printgml

    def use_rom(self):
        self.use_rom = True

    def set_work_dir(self):
        if not self.use_rom:
            self.work_dir = "/dev/shm/" + self.work_dir

    def get_max_num_repr(self):
        return self.max_num_repr

    def set_max_num_repr(self, max_num_repr):
        self.max_num_repr = max_num_repr

    def get_blast_thr(self):
        return self.blast_thr

    def set_blast_thr(self, blast_thr):
        self.blast_thr = blast_thr

    def get_threads(self):
        return self.threads

    def set_threads(self, threads):
        self.threads = threads

    def get_dom_file(self):
        return self.dom_file

    def set_dom_file(self, dom_file):
        self.dom_file = dom_file

    def get_blast_file(self):
        return self.blast_file

    def set_blast_file(self, blast_file):
        self.blast_file = blast_file

    def get_out(self):
        return self.out

    def set_out(self, out):
        self.out = out

    def get_num_threads_for_external_process(self):
        return self.num_threads_for_external_process

    def set_num_threads_for_external_process(self, available_cores):
        if self.threads <= 2 * available_cores / 3:
            self.num_threads_for_external_process = (2 * available_cores) // (3 * self.threads)


class SimGIC:
    def __init__(self):
        self.graph = None

    def set_graph_owl(self, graph):
        self.graph = graph

    def compute_distance(self, list1, list2, simlimit=None):
        parentset1 = set()
        parentset2 = set()

        if isinstance(list1, list):
            for n in list1:
                parentset1.update(self.graph.get_all_go_ancestors(n))
        else:
            parentset1.update(self.graph.get_all_go_ancestors(list1))

        if isinstance(list2, list):
            for n in list2:
                parentset2.update(self.graph.get_all_go_ancestors(n))
        else:
            parentset2.update(self.graph.get_all_go_ancestors(list2))

        dist = self.similarity(parentset1, parentset2)
        if simlimit is None:
            return dist
        else:
            if dist >= simlimit:
                return dist
            else:
                return 0.0

    @staticmethod
    def get_max_distance():
        return 1.0

    @staticmethod
    def get_min_distance():
        return 0.0

    @staticmethod
    def similarity(parentset1, parentset2):
        union = set()
        sumunion = 0.0
        suminter = 0.0

        union.update(parentset1)
        union.update(parentset2)

        intersection = parentset1.intersection(parentset2)

        if len(intersection) > 0:
            for n in intersection:
                suminter += n.getIC()

            for n in union:
                sumunion += n.getIC()

            # Not true if and only if parents1 and parents2 contain only the root
            if sumunion > 0:
                return suminter / sumunion

        return 0.0


class SSNetwork:
    def __init__(self, query, sg):
        self.query = query
        self.modularity = None

        insp = ConnectivityInspector(sg)
        vertexconnected = insp.connected_sets()
        for vertices in vertexconnected:
            subg = AsSubgraph(sg, vertices)
            map = {e: sg.get_edge_weight(e) for e in subg.edge_set()}
            tG = AsWeightedGraph(GraphTypeBuilder().weighted(True).undirected().build(), subg, map)
            if tG.contains_vertex(query):
                self.modularity = Modularity(tG, str)
                break
        else:
            self.modularity = None

    def get_all_clusters(self):
        clusters = []
        if self.modularity:
            for g in self.modularity.get_clusters():
                clusters.append(g.vertex_set())
        return clusters

    def get_cluster(self):
        cluster = set()
        if self.modularity:
            for g in self.modularity.get_clusters():
                if g.contains_vertex(self.query):
                    cluster = g.vertex_set()
                    break
        return cluster


class ThreadData:
    def __init__(self, prot_custers, neighbor_map):
        self.prot_custers = prot_custers
        self.neighbor_map = neighbor_map

    def get_prot_custers(self):
        return self.prot_custers

    def set_prot_custers(self, prot_custers):
        self.prot_custers = prot_custers

    def get_neighbor_map(self):
        return self.neighbor_map

    def set_neighbor_map(self, neighbor_map):
        self.neighbor_map = neighbor_map


class Tuple:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __eq__(self, other):
        if other is None:
            return False
        if other is self:
            return True
        if not isinstance(other, Tuple):
            return False
        return other.x == self.x and other.y == self.y

    def __hash__(self):
        prime = 31
        result = 1
        result = prime * result + (0 if self.x is None else hash(self.x))
        result = prime * result + (0 if self.y is None else hash(self.y))
        return result


if __name__ == '__main__':
    domains = Domains()
    domains.run_domains()

    delete_directory(domains.get_settings().get_work_dir())
