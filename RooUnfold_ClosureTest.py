
import sys
method = "bayes"
if len(sys.argv) > 1: method = sys.argv[1]

import ROOT
from ROOT import gRandom, TH1, TH1D, TCanvas
import atlasplots as aplt
import math
from ROOT import RooUnfoldResponse
try:
  cout= ROOT.std.cout  
except:
  cout= ROOT.cout

#  open MC file and get hTrue  hMeas and mig matrix

MC_Input1   = ROOT.TFile.Open("MCdis.root", "READ");
hTrue= MC_Input1.Get("TruthSelection/LpT_Truth_cut4");
hMeas= MC_Input1.Get("WplusenuSelection/elPt_cut7");
mig_matrix = MC_Input1.Get("WplusenuSelection/LpT_Reco_v_Truth_cut7");

# get the Acceptance and Efficiency 

Measured_Projec = mig_matrix.ProjectionX()
Truth_Projec    = mig_matrix.ProjectionY()

Efficiency_hist = Measured_Projec.Clone("Acceptance_hist")
Acceptance_hist = Truth_Projec.Clone("Efficiency_hist")

Acceptance_hist.Divide(hTrue)
Efficiency_hist.Divide(hMeas)

# Measured corrected
for i in range(0, hMeas.GetNbinsX()):
  hMeas.SetBinContent(i+1, hMeas.GetBinContent(i+1)*Efficiency_hist.GetBinContent(i+1))

ROOT.gROOT.SetBatch(True)

# Response matrix creation
response = RooUnfoldResponse(0,0,mig_matrix,"UNFOLD","UNFOLD");

# print resp matrix
response.Print()

# Plot response matrix
xmin = 25
xmax = 80
ymin = 25
ymax = 100
R = response.HresponseNoOverflow();
c7 = ROOT.TCanvas();

R.GetXaxis().SetRangeUser(xmin, xmax)
R.GetYaxis().SetRangeUser(ymin, ymax)

R.SetStats(0);
R.Draw("colz");
c7.Draw();
c7.SaveAs("response-matrix.png");

# Bayesian Unfolding
if   method == "bayes":
  unfold= ROOT.RooUnfoldBayes (response, hMeas, 4);#  OR             

hReco= unfold.Hreco();        

unfold.PrintTable (cout, hTrue);

#hReco correction
for i in range(0, hReco.GetNbinsX()):

    bin_content = Acceptance_hist.GetBinContent(i+1)
    if bin_content != 0:
        hReco.SetBinContent(i+1, hReco.GetBinContent(i+1) / bin_content)
    else:    
        hReco.SetBinContent(i+1, 0)  

#get the covriance matrix
covMatrix = unfold.Ereco(2);
print("Dimension de la matrice de covariance bayes:", covMatrix.GetNrows(), "x", covMatrix.GetNcols());
print("Matrice de covariance Bayes :");
covMatrix.Print()

# Set the ATLAS Style
aplt.set_atlas_style()

# plot cov matrix-Create a figure and axes
fig, ax = aplt.subplots(1, 1, name="fig1", figsize=(800, 800))

num_bins_x = covMatrix.GetNcols()
num_bins_y = covMatrix.GetNrows()
hist = ROOT.TH2D("hist", "Covariance Matrix",num_bins_x, 0, 100, num_bins_y, 0, 100)

# fill the his with the covMtrix values
for i in range(num_bins_x):
    for j in range(num_bins_y):
        hist.SetBinContent(i+1, j+1, covMatrix[i][j])
hist.GetXaxis().SetRangeUser(25, 60)
hist.GetYaxis().SetRangeUser(25, 60)
ax.plot2d(hist,"COLZ"); 

# Change pad margins to allow space for z-axis colour bar and for ATLAS label
ax.set_pad_margins(right=0.20, top=0.08)

# Set axis titles
ax.set_xlabel("p_{T}^{e} [GeV]")
ax.set_ylabel("p_{T}^{e} [GeV]")

# Add the ATLAS Label
aplt.atlas_label(ax.pad.GetLeftMargin(), 0.97, text="Internal", align=13)
ax.text(1 - ax.pad.GetRightMargin(), 0.97, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=33)

fig.savefig("covariance matrix.png")

# correlation matrix calculation
correlationMatrix = covMatrix.Clone()
for i in range(correlationMatrix.GetNrows()):
    for j in range(correlationMatrix.GetNcols()):
        if covMatrix[i][i] != 0.0 and covMatrix[j][j] != 0.0:
            correlationMatrix[i][j] = covMatrix[i][j] / (math.sqrt(covMatrix[i][i]) * math.sqrt(covMatrix[j][j]))
        else:
            correlationMatrix[i][j] = 0.0
correlationMatrix.Print()

#plot corr matrix- Create a figure and axes
fig1, ax1 = aplt.subplots(1, 1, name="fig2", figsize=(800, 800))
numb_bins_x = correlationMatrix.GetNcols()
numb_bins_y = correlationMatrix.GetNrows()
hist1 = ROOT.TH2D("hist", "Covariance Matrix",numb_bins_x, 0, 100, numb_bins_y, 0, 100)

# Fill the hist with correlation mateix values
for i in range(numb_bins_x):
    for j in range(numb_bins_y):
        hist1.SetBinContent(i+1, j+1, correlationMatrix[i][j])

hist1.GetXaxis().SetRangeUser(25, 60)#serangeruser --zoom
hist1.GetYaxis().SetRangeUser(25, 60)

ax1.plot2d(hist1, "COLZ")
# Change pad margins to allow space for z-axis colour bar and for ATLAS label
ax1.set_pad_margins(right=0.20, top=0.08)

# Set axis titles
ax1.set_xlabel("p_{T}^{e} [GeV]")
ax1.set_ylabel("p_{T}^{e} [GeV]")
#ax.set_zlabel("Events / (0.2 #times 5)", titleoffset=1.2)

# Add the ATLAS Label
aplt.atlas_label(ax1.pad.GetLeftMargin(), 0.97, text="Internal", align=13)
ax.text(1 - ax.pad.GetRightMargin(), 0.97, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=33)
fig1.savefig("correlation matrix.png")

# x-axis adjustment hMeas 
conversion_factor = 0.001

# Convert hMeas x-axis limits from MeV to GeV
hMeas.GetXaxis().SetLimits(hMeas.GetXaxis().GetXmin() * conversion_factor, hMeas.GetXaxis().GetXmax() * conversion_factor)

hRecoGraph = aplt.root_helpers.hist_to_graph(hReco)
hRecoGraph.SetLineWidth(0)
hRecoGraph.SetMarkerStyle(21)
hRecoGraph.SetMarkerSize(0.9)
hRecoGraph.SetMarkerColor(7)

# plot hTrue hMeas hReco with the ratio
fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(930, 930), hspace=0.05)

# Draw the histograms on these axes

ax1.plot(hRecoGraph,"LEP",label="Unfolded distribution", labelfmt="EP X0")
hMeas.SetLineStyle(5)
ax1.plot(hMeas, LineColor=3, label="Measured distribution", labelfmt="L")
hTrue.SetLineStyle(2)
ax1.plot(hTrue, linecolor=2, label="True distribution", labelfmt="L")
ax1.set_ylim(0, 40000)
ax1.set_xlim(25, 65) 

# Draw line at y=1 in ratio panel
line = ROOT.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
ax2.plot(line)

# Calculate and draw the ratio
ratio_hist = hReco.Clone("ratio_hist")
ratio_hist1 = hTrue.Clone("ratio_hist")
ratio_hist.Divide(hTrue)
ratio_hist1.Divide(hTrue)
ax2.plot(ratio_hist, linecolor=7,linewidth=2)#linecolor=ROOT.kRed+1
ax2.plot(ratio_hist1, linecolor=2)
# Add extra space at top of plot to make room for labels
ax1.add_margins(top=0.21)

# Set axis titles
ax1.set_ylabel("Events ", titleoffset=2.2)
ax2.set_ylabel("Unfold/True", loc="centre")
ax2.set_xlabel("p_{T}^{e} [GeV]", titleoffset=1)
ax2.set_xlim(25, 65)
ax2.set_ylim(0.5,1.5)  

# Go back to top axes to add labels
ax1.cd()

# Add the ATLAS Label
aplt.atlas_label(text="Internal", loc="upper left")
ax1.text(0.2, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=13)

# Add legend
ax1.legend(loc=(0.78, 0.74, 0.6, 0.90))

# Save the plot as a Png
fig.savefig("RooUnfold-pTe.png")

''''
##############PLOT EFFECIENCY AND ACCEPATNCE HIST#################
# Create a figure and axes
fig7, ax = aplt.subplots(1, 1, name="fig1", figsize=(800, 600))
# Draw the histograms on the same axis
ax.plot(Efficiency_hist, linecolor=ROOT.kOrange+1, label="Efficiency", linewidth=2, labelfmt="L", size=18)
ax.plot(Acceptance_hist, linecolor=ROOT.kAzure+1, label="Acceptance",linewidth=2, labelfmt="L", size=18)
ax.set_xlim(25,60)
ax.set_ylim(0, 1.2)

# Add extra space at top of plot to make room for labels
ax.add_margins(top=0.16)

# Set axis titles
ax.set_ylabel("Correction factors", titleoffset=2.2,size=18)
ax.set_xlabel("p_{T}^{e}[GEV]", titleoffset=1,size=18)

# Go back to top axes to add labels
ax.cd()

# Add the ATLAS Label
aplt.atlas_label(text="Internal", loc="upper left")
ax.text(0.2, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=13)

# Add legend
ax.legend(loc=(0.78, 0.78, 0.6, 0.90))

# Save the plot as a Png
fig7.savefig("Efficiency-acceptance.png")

