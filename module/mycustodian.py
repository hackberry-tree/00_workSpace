#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
自分用に修正
"""
# pylint: disable=unexpected-keyword-arg, no-value-for-parameter
from custodian.vasp.handlers import VaspErrorHandler
from custodian.utils import backup
from pymatgen.io.vaspio.vasp_input import VaspInput
from pymatgen.transformations.standard_transformations \
    import SupercellTransformation

from custodian.ansible.interpreter import Modder
from custodian.ansible.actions import FileActions, DictActions

class myVaspErrorHandler(VaspErrorHandler):
    """
    VaspErrorHandlerの修正版
    """
    def __init__(self, output_file="vasp.out"):
        VaspErrorHandler.__init__(self, output_file)

    def correct(self): #pylint: disable=R0915,R0912,R0914
        """
        "brmix"にてINCARにADDGRIDを追加する
        """
        backup([self.output_filename, "INCAR", "KPOINTS", "POSCAR", "OUTCAR",
                "vasprun.xml"])
        actions = []
        vi = VaspInput.from_directory(".")

        if "tet" in self.errors or "dentet" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ISMEAR": 0}}})

        if "inv_rot_mat" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"SYMPREC": 1e-8}}})

        if "brmix" in self.errors and "NELECT" in vi["INCAR"]:
            #brmix error always shows up after DAV steps if NELECT is specified
            self.errors.remove("brmix")
        if "brmix" in self.errors or "zpotrf" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ISYM": 0}}})
            #140814 add
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"ADDGRID": True}}})
            # Based on VASP forum's recommendation, you should delete the
            # CHGCAR and WAVECAR when dealing with these errors.
            actions.append({"file": "CHGCAR",
                            "action": {"_file_delete": {'mode': "actual"}}})
            actions.append({"file": "WAVECAR",
                            "action": {"_file_delete": {'mode': "actual"}}})

        if "subspacematrix" in self.errors or "rspher" in self.errors or \
                        "real_optlay" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"LREAL": False}}})

        if "tetirr" in self.errors or "incorrect_shift" in self.errors or \
                    "rot_matrix" in self.errors:
            actions.append({"dict": "KPOINTS",
                            "action": {"_set": {"generation_style": "Gamma"}}})

        if "amin" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"AMIN": "0.01"}}})

        if "triple_product" in self.errors:
            s = vi["POSCAR"].structure
            trans = SupercellTransformation(((1, 0, 0), (0, 0, 1), (0, 1, 0)))
            new_s = trans.apply_transformation(s)
            actions.append({"dict": "POSCAR",
                            "action": {"_set": {"structure": new_s.to_dict}},
                            "transformation": trans.to_dict})

        if "pricel" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"SYMPREC": 1e-8, "ISYM": 0}}})

        if "brions" in self.errors:
            potim = float(vi["INCAR"].get("POTIM", 0.5)) + 0.1
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"POTIM": potim}}})

        if "zbrent" in self.errors:
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"IBRION": 1}}})

        if "too_few_bands" in self.errors:
            if "NBANDS" in vi["INCAR"]:
                nbands = int(vi["INCAR"]["NBANDS"])
            else:
                with open("OUTCAR") as f:
                    for line in f:
                        if "NBANDS" in line:
                            try:
                                d = line.split("=")
                                nbands = int(d[-1].strip())
                                break
                            except (IndexError, ValueError):
                                pass
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"NBANDS": int(1.1 * nbands)}}})

        if "aliasing" in self.errors:
            with open("OUTCAR") as f:
                grid_adjusted = False
                changes_dict = {}
                for line in f:
                    if "aliasing errors" in line:
                        try:
                            grid_vector = line.split(" NG", 1)[1]
                            value = [int(s) for s in grid_vector.split(" ")
                                     if s.isdigit()][0]

                            changes_dict["NG" + grid_vector[0]] = value
                            grid_adjusted = True
                        except (IndexError, ValueError):
                            pass
                    #Ensure that all NGX, NGY, NGZ have been checked
                    if grid_adjusted and 'NGZ' in line:
                        actions.append({"dict": "INCAR",
                                        "action": {"_set": changes_dict}})
                        break

        m = Modder(actions=[DictActions, FileActions])
        modified = []
        for a in actions:
            if "dict" in a:
                modified.append(a["dict"])
                vi[a["dict"]] = m.modify_object(a["action"], vi[a["dict"]])
            elif "file" in a:
                m.modify(a["action"], a["file"])
        for f in modified:
            vi[f].write_file(f)
        return {"errors": list(self.errors), "actions": actions}
