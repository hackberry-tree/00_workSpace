#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
list_run.txtに記載のjobを逐次runするscript
"""
import os
import re
import sys
import getpass
from subprocess import Popen, PIPE
from commopy import Cabinet

def main():
    """
    引数でjobの数を指定する
    """
    try:
        num = int(sys.argv[1])
    except IndexError:
        num = 1
    job = JobHandler(num_run=num)
    job.exe_batch()

class JobHandler(object):
    """
    jobを管理する
    run_fileには、"cd $PBS_O_WORKDIR"をpathに置き換える処理を行うので
    "cd $PBS_O_WORKDIR"の記述が必須
    """
    def __init__(self, num_run=1, fname_run='list_run',
                 fname_running='running_jobs'):
        self.num_run = num_run
        self.list_run_file = fname_run
        self.running_jobs_file = fname_running

        lines = Cabinet.read_file(fname_run)
        self.finished_list = [x.split() for x in lines if x[0] == '#']
        self.run_list = [x.split() for x in lines if x[0] != '#']
        try:
            lines = Cabinet.read_file(fname_running)
            self.running_jobs = [x for x in lines if x != '\n']
        except IOError:
            self.running_jobs = []

    def exe_batch(self):
        """
        batch jobを実行する
        """
        self.is_running()
        while self.num_run > self.num_running:
             if not self.run_list:
                print('false')
                return
             self.prep_run_file()
             self.run_job_and_get_id()
             self.is_running()
        print('true')


    @property
    def num_running(self):
        """
        動作中のjobの数を出力
        """
        return len(self.running_jobs)

    def prep_run_file(self, next_run='next_run.sh'):
        """
        run_fileを修正してnext_run.shを作成
        run_listの実行したjobをコメントアウト
        """
        work_path, fname_exe = self.run_list[0]
        key = "cd $PBS_O_WORKDIR\n"
        # work_pathが相対pathの場合、絶対pathに修正する
        if work_path[0] == '/':
            alt = "cd {0}\n".format(os.path.join(work_path))
        else:
            alt = "cd {0}\n".format(os.path.join(os.getcwd(), work_path))

        lines = Cabinet.read_file(fname_exe)

        # $PBS_O_WORKDIRの記述が無い場合errorメッセージを出力
        try:
            pos = lines.index(key)
        except ValueError:
            print("{0}に'cd $PBS_O_WORKDIR'の記述がありません".format(fname_exe))
            exit()
        lines[pos] = alt
        Cabinet.write_file(next_run, lines)

        finished = self.run_list.pop(0)
        finished[0] = '#' + finished[0]
        self.finished_list.append(finished)
        tmp_list = self.finished_list + self.run_list
        all_list = [" ".join(x) + "\n" for x in tmp_list]
        Cabinet.write_file(self.list_run_file, all_list)

    def run_job_and_get_id(self, fname_exe='next_run.sh'):
        """
        jobを走らせる
        走らせるファイルはnext_run.sh
        """
        cmd = ['qsub', fname_exe]
        job_id = Popen(cmd, stdout=PIPE).communicate()[0]
        add_line = job_id.split('.')[0] + "\n"
        self.running_jobs.append(add_line)
        Cabinet.write_file(self.running_jobs_file, self.running_jobs)

    def is_running(self):
        """
        jobが動作中かどうかcheckする
        終了していた場合、running_jobsから削除
        """
        cmd = ['qstat']
        jobs = Popen(cmd, stdout=PIPE).communicate()[0].split('\n')
        uname = getpass.getuser()
        still_run = []
        for job_id in self.running_jobs:
            key = r"{0}.*{1}.*".format(job_id[:-1], uname)
            meta = re.compile(key)
            judge = True
            for job in jobs:
                if meta.match(job):
                    judge = None
            if not judge:
                still_run.append(job_id)
        self.running_jobs = still_run
        Cabinet.write_file(self.running_jobs_file, self.running_jobs)


if __name__ == "__main__":
    main()
