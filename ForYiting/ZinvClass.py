import ROOT
import argparse
import math

from RA1Object import RA1Object
from collections import OrderedDict
from tdrStyle import SetPalette

parser = argparse.ArgumentParser(description='parser for Zinvclass')
parser.add_argument('--inputDir',action="store",help="input directory")
parser.add_argument('--inputPath',action="store",help="input path")
parser.add_argument('--outputDir',action="store",help="output directory")
parser.add_argument('--outputPath',action="store",help="output path")
parser.add_argument('--verbose',action="store_true")
parser.add_argument("--what",action="store",default="fit_b")
parser.add_argument("--noWrite",action="store_true")

class FitReader(RA1Object):
    """
    Read results from Maxlikelihood fit, and register in analysis bin class
    """ 

    def __init__(self,name,**kargs):
        super(FitReader,self).__init__(name,**kargs)
        self.ctSystNames = ["alphaT","eq0b_eq1b_muon","muplus_muminus","phot_mumu","mu_mumu"]

    def readFitHist(self,analysisBins,mlfitFilePath,whichFit,region,histName,customName=None,verbose=False):
        self.mlfitFile = ROOT.TFile(mlfitFilePath,"READ")
        for anaBin in analysisBins:
            anaBin.mlFitFolderName = self.convertAnaBinStr(anaBin)
            attrName = customName if customName else whichFit+"_"+histName
            histPath = "shapes_"+whichFit+"/"+anaBin.mlFitFolderName+"_"+region+"/"+histName
            tempObj = self.mlfitFile.Get(histPath)
            if tempObj != None:
                setattr(anaBin,attrName,tempObj)
            else:
                setattr(anaBin,attrName,None)

    def getFitParams(self,analysisBins,mlfitFilePath,whichFit,verbose=False):
        """
        Read nuisance form Maxlikelihood fit, and register in analysis bin class
        """
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.mlfitFile = ROOT.TFile(mlfitFilePath,"READ")
        fitResult = self.mlfitFile.Get(whichFit)
        fitParams = fitResult.floatParsFinal()
        self.rateParamList = []
        self.correlatedSystList = []
        for i in range(fitParams.getSize()):
            nuis    = fitParams.at(i)
            name    = nuis.GetName()
            value   = nuis.getVal()
            error   = nuis.getError()
            if verbose:
                print name,value,error
            # If this nuis is rate parameter
            if "Rate" in name:
                rateParamName = self.registerRateParams(analysisBins,name,value,error)
                self.rateParamList.append(rateParamName)
            # if this nuis is a CT nuis
            elif self.isCTNuis(name):
                self.registerCTNuis(analysisBins,name,value,error,mergeHtBin="mumu" in name)
            elif "qcdSyst" in name:
                pass
            else:
                for anaBin in analysisBins:
                    setattr(anaBin,name,value)
                    setattr(anaBin,name+"Err",error)
                if name not in self.correlatedSystList:
                    self.correlatedSystList.append(name)
            # else:
            #     print name
        self.registerEmptyBins(analysisBins,self.rateParamList)                
        self.registerEmptyBins(analysisBins,self.ctSystNames)     
        self.registerEmptyBins(analysisBins,self.correlatedSystList)     

    @staticmethod
    def registerCTNuis(analysisBins,name,value,error,mergeHtBin=False):
        for anaBin in analysisBins:
            if not mergeHtBin:
                anaBinStr = "_".join([anaBin.htBinStr(True),anaBin.nJetCatStr()])
            else:
                anaBinStr = "_".join([anaBin.htBinStr(capitalize=True,mergeHtBin=True),anaBin.nJetCatStr()])
            if anaBinStr in name:
                attrName = name.replace("_"+anaBinStr,"")
                setattr(anaBin,attrName,value)
                setattr(anaBin,attrName+"Err",error)

    def isCTNuis(self,name):
        for systName in self.ctSystNames:
            if systName in name:
                return systName
        return None

    @staticmethod
    def registerEmptyBins(analysisBins,paramList):
        for paramName in paramList:
            for anaBin in analysisBins:
                if not hasattr(anaBin,paramName):
                    setattr(anaBin,paramName,0.)
                    setattr(anaBin,paramName+"Err",0.)

    @staticmethod
    def registerRateParams(analysisBins,name,value,error):
        for anaBin in analysisBins:
            anaBinStr = "_".join([anaBin.bJetStr,anaBin.nJetStr,anaBin.htBinStr(True)])
            if anaBinStr in name:
                #print anaBin,name
                #print value
                #print "------------------"
                attrName = name.replace(anaBinStr,"")
                setattr(anaBin,attrName,value)
                setattr(anaBin,attrName+"Err",error)
                break
        return attrName

    @staticmethod
    def convertAnaBinStr(anaBin):
        """Convert str(anaBin) to htBin_nBJetBin_nJet structure, the mlfit folder name format"""
        nJetStr = anaBin.nJetStr if anaBin.nJetStr != None else "Inc"
        bJetStr = anaBin.bJetStr if anaBin.bJetStr != None else "Inc"
        return "ht"+anaBin.htBinStr(True)+"_"+bJetStr+"_"+nJetStr

    @staticmethod
    def convertHtCatDict(htCatDict):
        tempHtBinDict = OrderedDict()
        for cat,htBins in htCatDict.iteritems():
            bjet = cat.bjet if cat.bjet != "Inc" else None
            njet = cat.njet if cat.njet != "Inc" else None
            tempHtBinDict[(bjet,njet)] = htBins
        return tempHtBinDict

class FitWrapper(RA1Object):
    @staticmethod
    def makeSummary(analysisBins,attrName=None,lambdaContentFunc=None,lambdaErrorFunc=None,histName="summary"):
        """
        Make summary plot in (nJet,nB) vs HT format, based on analysis bin attribute, determined by attrName and lambda string
        Aim to standardize the process to make summary
        """
        maxHtBins = max([len(anaBins) for jetCat,anaBins in analysisBins.binDict.iteritems() ])
        summaryHist = ROOT.TH2D(histName,histName,maxHtBins,-0.5,maxHtBins-0.5,len(analysisBins.binDict),-0.5,len(analysisBins.binDict)-0.5)
        iJetCat = 0
        for jetCat,anaBins in analysisBins.binDict.iteritems():
            iJetCat += 1
            nJetStr = jetCat[1] if jetCat[1] != None else "Inc"
            nBJetStr = jetCat[0] if jetCat[0] != None else "Inc"
            jetCatStr = nBJetStr+"_"+nJetStr
            summaryHist.GetYaxis().SetBinLabel(iJetCat,jetCatStr)
            for iHt,anaBin in enumerate(anaBins):
                summaryHist.GetXaxis().SetBinLabel(iHt+1,str(anaBin.htBin[0]))
                if lambdaContentFunc and lambdaErrorFunc:
                    binContent = lambdaContentFunc(anaBin)
                    binError = lambdaErrorFunc(anaBin)
                else:
                    binContent = getattr(anaBin,attrName)
                    binError = getattr(anaBin,attrName+"Err")
                summaryHist.SetBinContent(iHt+1,iJetCat,binContent)
                summaryHist.SetBinError(iHt+1,iJetCat,binError)

        return summaryHist  

    @staticmethod
    def emptyHist(hist):
        nBinsX = hist.GetXaxis().GetNbins()
        nBinsY = hist.GetYaxis().GetNbins()
        for xbin in range(1,nBinsX+1):
            for ybin in range(1,nBinsY+1):
                hist.SetBinContent(xbin,ybin,0.)
                hist.SetBinError(xbin,ybin,0.)

    @staticmethod
    def make2DPull(hist1,hist2,pullHistName="",combineError=False,verbose=False,pullCut=3.):
        """
        make 2D pull distribution based on two hists with the same format
        error will be taken from hist2: (content2-content1)/error2
        """
        nBinsX = hist1.GetXaxis().GetNbins()
        nBinsY = hist1.GetYaxis().GetNbins()
        nBinsMiss = 0
        if pullHistName:
            pullHist = hist1.Clone(pullHistName)
        else:
            pullHist = hist1.Clone("pull_"+hist1.GetName())
        FitWrapper.emptyHist(pullHist)
        for xbin in range(1,nBinsX+1):
            for ybin in range(1,nBinsY+1):
                binContent1 = hist1.GetBinContent(xbin,ybin)
                binContent2 = hist2.GetBinContent(xbin,ybin)
                error = hist2.GetBinError(xbin,ybin)**2
                xLabel = hist2.GetXaxis().GetBinLabel(xbin)
                yLabel = hist2.GetYaxis().GetBinLabel(ybin)
                if combineError:
                    error += hist1.GetBinError(xbin,ybin)**2
                if binContent1 and binContent2 and error:
                    pull = (binContent2-binContent1)/math.sqrt(error)
                    pullHist.SetBinContent(xbin,ybin,pull)
                    pullHist.SetBinError(xbin,ybin,0.)
                    if abs(pull) > pullCut and verbose:
                        print "bin with large pull (> %s): "%pullCut, xLabel, yLabel
                elif binContent1 and not binContent2:
                    nBinsMiss += 1
                    if verbose:
                        print "Missing hist2 bins", xLabel, yLabel
        print "Number of missing bins: ", nBinsMiss
        return pullHist

    @staticmethod
    def make1DPull(hist1,hist2,pullHistName="",combineError=False,verbose=False):
        """
        make 1D pull distribution based on two hists with the same format
        error will be taken from hist2: (content2-content1)/error2
        """
        nBinsX = hist1.GetXaxis().GetNbins()
        nBinsY = hist1.GetYaxis().GetNbins()
        pull1DHist = ROOT.TH1D(pullHistName,pullHistName,44,-10.5,10.5)
        for xbin in range(1,nBinsX+1):
            for ybin in range(1,nBinsY+1):
                binLabelY = hist1.GetYaxis().GetBinLabel(ybin)
                #if "eq0b" in binLabelY or "eq1b" in binLabelY: continue
                #if "eq2b" in binLabelY or "ge3b" in binLabelY: continue
                binContent1 = hist1.GetBinContent(xbin,ybin)
                binContent2 = hist2.GetBinContent(xbin,ybin)
                error = hist2.GetBinError(xbin,ybin)**2
                if combineError:
                    error += hist1.GetBinError(xbin,ybin)**2
                if binContent1 and binContent2 and error:
                    pull1DHist.Fill((binContent2-binContent1)/math.sqrt(error))
        return pull1DHist

    
class DrawSetting(object):
    def __init__(self):
        self.generalSetting = {"textFormat":"4.2f","Stats":"nemr",}
        self.th1Setting = {"Stats":0,}
        self.th2Setting = {"Stats":0,}

defaultDrawSetting = DrawSetting()

class Drawer(object):
    @staticmethod
    def drawEverything(histDict,outputDir,setting=defaultDrawSetting):
        setPalette = SetPalette()
        setPalette('kBird')
        ROOT.gStyle.SetPaintTextFormat(setting.generalSetting["textFormat"])
        ROOT.gStyle.SetOptStat(setting.generalSetting["Stats"])
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        c = ROOT.TCanvas()
        for fileName,itemList in histDict.iteritems():
            c.Print(outputDir+"/"+fileName+"[")
            for item in itemList:
                className = item.__class__.__name__
                if className == "tuple" or className == "list":
                    maximum = max( [ hist.GetMaximum() for hist in item ] )
                    minimum = max( [ hist.GetMinimum() for hist in item ] )
                    leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
                    for i,hist in enumerate(item):
                        hist.GetYaxis().SetRangeUser(0.8*minimum,1.2*maximum)
                        hist.SetStats(0)
                        hist.SetTitle(hist.GetName())
                        hist.SetLineColor(i+2)
                        leg.AddEntry(hist,hist.GetName())
                        if i == 0:
                            hist.DrawCopy()
                        else:
                            hist.DrawCopy("same")
                    leg.Draw()
                elif "TH1" in className:
                    item.SetStats(setting.th1Setting["Stats"])
                    if "setYRange" in setting.th1Setting:
                        item.GetYaxis().SetRangeUser(*setting.th1Setting["setYRange"])
                    if "markerSize" in setting.th1Setting:
                        item.SetMarkerSize(setting.th1Setting["markerSize"])

                    item.SetTitle(item.GetName())
                    if "drawOption" in setting.th1Setting:
                        item.DrawCopy(setting.th1Setting["drawOption"])
                    else:
                        item.DrawCopy("TEXT E1")
                elif "TH2" in className:
                    item.SetStats(setting.th2Setting["Stats"])
                    if "setZRange" in setting.th2Setting:
                        item.GetZaxis().SetRangeUser(*setting.th2Setting["setZRange"])
                    if "markerSize" in setting.th2Setting:
                        item.SetMarkerSize(setting.th2Setting["markerSize"])
                    item.SetTitle(item.GetName())
                    if "drawOption" in setting.th2Setting:
                        item.DrawCopy(setting.th2Setting["drawOption"])
                    else:
                        item.DrawCopy("colztexte")
                c.Print(outputDir+"/"+fileName)
            c.Print(outputDir+"/"+fileName+"]")

