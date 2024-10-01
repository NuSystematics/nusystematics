import ROOT

f_out = ROOT.TFile("ZExpReweighTemplate.root", "RECREATE")

h_rw = ROOT.TH2D("h_rw_Enu_vs_Q2", "", 10, 0., 10., 10, 0., 2.)

for ix in range(0, 10):
  ibinx = ix+1
  for iy in range(0, 10):
    ibiny = iy+1

    h_rw.SetBinContent(ibinx, ibiny, 1.)

f_out.cd()
h_rw.Write()
f_out.Close()
