#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
test = "tot        -0.03   -0.07    7.48    7.39"

keywords = r"\s*([\d]+)\s*F=\s*([\d.]+E\+[\d.]+)\s*E0=\s*.*"
keywords = r"\s*tot\s+[\d\-\.]+\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)\s*.*"
keywords = (r"\s*tot\s+[\d\-\.]+\s+[\d\-\.]+"
            r"\s+[\d\-\.]+\s+([\d\-\.]+)\s*.*")
meta = re.compile(keywords)


print(meta.match(test))
print(meta.match(test).group(1))
