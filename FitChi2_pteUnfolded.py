import ROOT 
import ROOT as root
import atlasplots as aplt 
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TF1, TGraph, gROOT, gRandom
import numpy as np
								
# Open the ROOT file

Unfolded_file = ROOT.TFile.Open("UNFOLDED_DISTRIBUTIONS.root")

# Open the MC nominal hists
McNom_hist = Unfolded_file.Get("hReco_MCnominal")
Cov_Nom = Unfolded_file.Get("CovMatrix_MCnominal")
	
# Define a list of unfolded hist (Nom+templates)
hist_names = ["hReco_MC-200","hReco_MC-150","hReco_MC-100","hReco_MC-50","hReco_MC-25","hReco_MC0","hReco_MC25","hReco_MC50","hReco_MC100","hReco_MC150","hReco_MC200"]

# Define a list of covriance matrices
cov_names = ["CovMatrix_MC-200","CovMatrix_MC-150","CovMatrix_MC-100","CovMatrix_MC-50","CovMatrix_MC-25","CovMatrix_MC0","CovMatrix_25","CovMatrix_50","CovMatrix_100","CovMatrix_150","CovMatrix_200"]

# Conversion of histogrammes to tables numpy:
# Mcnominal as array
McNom_array = np.empty((Cov_Nom.GetNbinsX()))
# Loop over histogram bins and assign values to NumPy array
for i in range(1, McNom_hist.GetNbinsX() + 1):
	McNom_array[i-1] = McNom_hist.GetBinContent(i)
McNom_array = McNom_array[25:]

for i in range(len(McNom_array)):
    contenu_bin = McNom_array[i]
    print("Contenu du bin", i+1, ":", contenu_bin)

# MC_nominal_covariance as array
MCnom_covariance_array = np.empty((Cov_Nom.GetNbinsX(), Cov_Nom.GetNbinsY()))

# Loop over histogram bins and assign values to NumPy array
for i in range(1, Cov_Nom.GetNbinsX() + 1):
	for j in range(1, Cov_Nom.GetNbinsY() + 1):
		MCnom_covariance_array[i-1, j-1] = Cov_Nom.GetBinContent(i, j)

MCnom_covariance_array = MCnom_covariance_array[25:, 25:] 

#MC templates as array
unfolded_arrays = []
unfolded_covariance_arrays = []

# Loop over the histograms
for i in range(len(hist_names)):
	Unfold_hist = Unfolded_file.Get(hist_names[i])
	Unfold_hist_array= np.empty((Unfold_hist.GetNbinsX()))
	# Loop over histogram bins and assign values to NumPy array
	for j in range(1, Unfold_hist.GetNbinsX() + 1):
		Unfold_hist_array[j-1] = Unfold_hist.GetBinContent(j)
	Unfold_hist_array = Unfold_hist_array[25:]
	print(Unfold_hist_array.shape)    
	print(Unfold_hist_array)
	unfolded_arrays.append(Unfold_hist_array)

for k in range(len(cov_names)):    
	# get the histogram
	cov_hist = Unfolded_file.Get(cov_names[k])
	unfolded_covariance_array = np.empty((cov_hist.GetNbinsX(), cov_hist.GetNbinsY()))
	# Loop over histogram bins and assign values to NumPy array
	for i in range(1, cov_hist.GetNbinsX() + 1):
		for j in range(1, cov_hist.GetNbinsY() + 1):
			unfolded_covariance_array[i-1, j-1] = cov_hist.GetBinContent(i, j)

	unfolded_covariance_array = unfolded_covariance_array[25:, 25:] 
	print(unfolded_covariance_array.shape) 

	unfolded_covariance_arrays.append(unfolded_covariance_array)
chi2list=[]

# chi2 calculation for each unfolded distribution
for i in range(len(unfolded_arrays)):
	# Calculating the difference between unfolded distributions
	diff = np.subtract(McNom_array, unfolded_arrays[i])
	print(diff)
	# Calculating the sum of covariance matrices
	total_covariance_array = np.add(MCnom_covariance_array,unfolded_covariance_arrays[i]) 
	print(total_covariance_array)
	# Chi2 calculation 
	chi2 = np.dot(np.dot(diff.T, np.linalg.inv(total_covariance_array)), diff)
	chi2list.append(chi2)

print(chi2list)

# Set the ATLAS Style
aplt.set_atlas_style()

# plot cov matrix-Create a figure and axes
fig, ax = aplt.subplots(1, 1, name="fig1", figsize=(800, 800))

#define mass step
MassSteps=[-0.2,-0.15,-0.1,-0.05,-0.025,0,0.025,0.05,0.1,0.15,0.2]

#define the graphs     
graph2 = TGraph(len(chi2list)) 
for i in range(0, len(chi2list)):
	for i in range(0, len(MassSteps)):
		graph2.SetPoint(i, 80.399+MassSteps[i], chi2list[i])  

# Fit function
fitchi22  = ROOT.TF1("fitchi22","(x-[0])*(x-[0])/[1]/[1] +[3]*x*x*x+ [4]*x*x*x*x +[2]",80.150,80.650) 
fitchi22.SetParameter(1,10)
fitchi22.SetParameter(2,5.)
fitchi22.SetParameter(3, 1.)
fitchi22.SetParameter(4, 1.)

#fit graph 2
graph2.Fit("fitchi22","rQ")
graph2.Fit("fitchi22","rQ")
graph2.Fit("fitchi22","rQ")

# print the results.
print("Nominal : ",   fitchi22.GetParameter(0))
print("Stat error : ",fitchi22.GetParameter(1))

results2 = []
results2.append(fitchi22.GetParameter(0))
results2.append(fitchi22.GetParameter(1))
print(results2)
c1 = TCanvas("c1","c1", 800, 800)
   
# Define x an y titles
graph2.GetXaxis().SetTitle("p_{T}^{e}[GeV]")  
graph2.GetYaxis().SetTitle("\chi^{2}")
graph2.GetYaxis().SetTitleOffset(1.2)
graph2.GetYaxis().SetTitleSize(0.05)
graph2.GetYaxis().SetTitleFont(42)

fitchi22.SetLineColor(ROOT.kRed)

graph2.SetLineWidth(0)
graph2.SetMarkerStyle(3)

# Access the y axis of the curve
yaxis = graph2.GetYaxis()

# Set y-axis range between 0 and 20
yaxis.SetRangeUser(0,400)

# draw graphs and fit functions in the same canvas
graph2.Draw()
fitchi22.Draw("same")
	
# Define left margin
c1.SetLeftMargin(0.15)

# Add titles
title1 = ROOT.TLatex()
title1.SetTextSize(0.030)
title1.SetTextAlign(11)
title1.SetTextFont(42)
title1.DrawLatexNDC(0.25, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}")
title2 = ROOT.TLatex()
title2.SetTextSize(0.036)
title2.SetTextAlign(11)
title2.SetTextFont(42)
title2.DrawLatexNDC(0.23, 0.88, "ATLAS Internal") 

# add a legend
legend = ROOT.TLegend(0.72, 0.82, 0.8, 0.91)
legend.AddEntry(graph2, "\chi^{2} values",)
legend.AddEntry(fitchi22, "Fit function", "l")
legend.SetTextSize(20)
legend.Draw()

#print c1
c1.SaveAs("chi2_Unfolded_MC.png")