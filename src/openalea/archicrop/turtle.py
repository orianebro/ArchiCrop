from __future__ import annotations

import openalea.plantgl.all as pgl


class CerealsTurtle(pgl.PglTurtle):
    def __init__(self):
        super().__init__()
        self.context = {}

    def transform(self, mesh, face_up=False):
        x = self.getUp()
        z = pgl.Vector3(0, 0, 1) if face_up else self.getHeading()
        bo = pgl.BaseOrientation(x, z ^ x)
        matrix = pgl.Transform4(bo.getMatrix())
        matrix.translate(self.getPosition())
        # print 'Position ', turtle.getPosition()
        return mesh.transform(matrix)

    def getFrame(self):
        pos = self.getPosition()
        Head = self.getHeading()
        Up = self.getUp()
        return {"Position": pos, "Head": Head, "Up": Up}

    def setFrame(self, frame):
        self.move(frame["Position"])
        self.setHead(frame["Head"], frame["Up"])
