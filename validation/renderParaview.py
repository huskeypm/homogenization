from paraview.simple import *
case = "0p50"
case = "1p50"
servermanager.LoadState(case+".pvsm")

try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

RenderView1 = GetRenderView()
WriteAnimation(case+'.png', Magnification=1, Quality=9, FrameRate=1.000000)
# doesn't work WriteImage('test.png',Quality=9)


