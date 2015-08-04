#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

with open("/Users/enoki/Desktop/stress.out", "r") as rfile:
    lines = rfile.readlines()

temp_dict = {x.split("/")[0]: [float(y) for y in x.split()[3:]] for x in lines}
str_name = []
stress = []

for key, st in temp_dict.items():
    str_name.append(key)
    stress.append(st)
str_name = np.array(str_name)
stress = np.array(stress)
print(len(stress))

mean = stress.mean(axis=-1) > 10
# print(mean.sum())

plus = (stress ** 2) ** 0.5
mean2 = plus.mean(axis=-1) > 35
print(mean2.sum())

mean3 = plus[:, 3:].mean(axis=-1) > 50
# print(mean3.sum())


with open("/Users/enoki/Desktop/energies.txt", "r") as rfile:
    lines_e = rfile.readlines()
out = ""
for l in lines_e:
    if l.split('*')[0] in str_name[mean2]:
        out += "*" + l
    else:
        out += l

print(out)

