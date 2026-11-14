# Simple portal frame analysed, saved, then read back: the forces read back are bit-for-bit identical.
import os
import tempfile
import keta

m = keta.Model()
a, b, c = m.add_node(0, 0, 0), m.add_node(0, 0, 4), m.add_node(6, 0, 4)
steel = m.add_material(210e9, 0.3, rho=7850.0)
section = m.add_section_rect(0.2, 0.4)
column, beam = m.add_beam(a, b, section, steel), m.add_beam(b, c, section, steel)
m.add_support(a, True, True, True, True, True, True)
m.add_support(c, False, True, True)
m.add_udl(beam, wz=-20e3, case=1)
m.add_point_load(b, fx=15e3, case=2)

r = m.analyze(stations=5)
path = os.path.join(tempfile.mkdtemp(), "frame.json")
r.save(path, note="portal frame, variant A")
info = keta.results_info(path)
print(f"schema {info['schema_version']}, note “{info['note']}”, {info['cases']} cases")

reread = m.read_results(path)
identical = all(reread.case(i).bar_forces(bar) == r.case(i).bar_forces(bar)
                for i in (1, 2) for bar in (column, beam))
print("forces read back identical:", identical)
print(f"moment at the column head, case 1: {reread.case(1).bar_forces(column)[-1][5] / 1e3:.2f} kN.m")
