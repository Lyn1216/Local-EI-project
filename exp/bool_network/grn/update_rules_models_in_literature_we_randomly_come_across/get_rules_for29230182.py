#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

a = '''ANG1(t + 1) =ANG1(t)
ANG2(t + 1) =(¬KLF2(t)) ∧ (HIF1(t) ∨ AP1(t) ∨ ETS(t) ∨ FOXO1(t)) 
TIE1(t + 1) =(Integrin(t) ∨ ANG1(t) ∨ ANG2(t))
TIE2(t + 1) =TIE1(t) ∧ ANG1(t) ∧ (¬ANG2(t))∧ ((¬VEPTP(t)) ∨ KLF2(t) ∨ ETS(t)) 
VEPTP(t + 1) =HIF2a(t) 
ShearStress(t + 1) =ShearStress(t)
mir126(t + 1) =KLF2(t) 
SPRED1(t + 1) =¬mir126(t) 
PI3KR2(t + 1) =¬mir126(t)
NF1(t + 1) =SPRED1(t)
KLF2(t + 1) =ShearStress(t) 
PIEZO1(t + 1) =ShearStress(t)
TRPV4(t + 1) =ShearStress(t) ∧ ¬cGMP(t) 
Integrin(t + 1) =(ShearStress(t) ∨ TIE2(t)) ∧ ETS(t)
GqG11(t + 1) =ShearStress(t) 
Actin(t + 1) =ShearStress(t)
V_Ecadherin(t + 1) =ETS(t) ∧ (SRC(t) ∨ FAK(t) ∨ Actin(t) ∨ VEPTP (t)) 
PECAM1(t + 1) =FYN(t) ∨ Actin(t) ∨ V_Ecadherin(t) 
SHP2(t + 1) =PECAM1(t)
FAK(t + 1) =SRC(t) ∨ Integrin(t)
SRC(t + 1) =FAK(t) ∨ TSAd(t) ∨ GqG11(t)
PI3K(t + 1) =AXL(t) ∨ (¬PI3KR2(t)) ∨ V_Ecadherin(t) ∨ TIE2(t)
AKT(t + 1) =PIP3(t)
TSC(t + 1) =AMPK(t) ∧ (¬AKT(t))
AXL(t + 1) =SRC(t)
PIP3(t + 1) =PI3K(t) ∧ (¬PTEN(t))
RHEB2(t + 1) =¬TSC(t)
FYN(t + 1) =TSAd(t)
TSAd(t + 1) =VEGFR22(t) ∨ VEGFR23(t)
PTEN(t + 1) =CSL(t) ∨ RHO(t)
BMP9(t + 1) =BMP9(t)
BMP10(t + 1) =BMP10(t)
TGFB1(t + 1) =TGFB1(t)
ACVR2A(t + 1) =1
BMPRII(t + 1) =1
TGFBRII(t + 1) =1
ALK1(t + 1) =(TGFB1(t) ∨ BMP9(t) ∨ BMP10(t))∧ (ACVR2A(t) ∨ BMPRII(t) ∨ TGFBRII(t))
ENG(t + 1) =TGFB1(t) ∨ BMP9(t) ∨ BMP10(t) 
ALK5(t + 1) =BMP9(t) ∧ TGFBRII(t)
BG(t + 1) =BMP9(t)
SMAD4(t + 1) =1
PDGFB(t + 1) =FOXO1(t) ∨ SMAD2(t)
TRAF6(t + 1) =ALK1(t)
TAK1(t + 1) =TRAF6(t)
sGC(t + 1) =1
cGMP(t + 1) =NO(t) ∧ sGC(t)
eNOS(t + 1) =AKT(t) ∨ SIRT1(t)
NO(t + 1) =eNOS(t) ∨ Calmodulin(t) 
Calcium(t + 1) =IP3(t) ∨ TRPV4(t) ∨ PIEZO1(t) 
Calmodulin(t + 1) =Calcium(t) 
Calcineurin(t + 1) =Calmodulin(t) 
NFAT(t + 1) =Calcineurin(t) 
ETS(t + 1) =ERK(t) 
RAS(t + 1) =SHP2(t) ∨ ALK1(t) ∨ (¬NF1(t))
RAF1(t + 1) =(PKCB2(t) ∨ RAS(t)) ∧ (¬AKT (t)) 
MEK(t + 1) =RAF1(t) ∨ FRS2a(t) 
ERK(t + 1) =MEK(t) ∨ VEGFR33(t) 
PLCg(t + 1) =VEGFR22(t) ∨ VEGFR33(t) ∨ DSH(t) ∨ FZD45(t)
IP3(t + 1) =PLCg(t) 
DAG(t + 1) =PLCg(t)
PKCB2(t + 1) =DAG(t) ∧ Calcium(t)
AMPATP(t + 1) =AMPATP(t) 
AMPK(t + 1) =(AMPATP(t) ∨ (¬Oxygen(t))) ∧ (¬AKT(t)) 
mTOR(t + 1) =AMPK(t) ∨ RHEB2(t) 
HIF1(t + 1) =mTOR(t) ∧ ¬PHDs(t) ∧ SIRT1(t) 
FOXO1(t + 1) =SIRT1(t) ∧ (¬AKT (t)) 
ART1(t + 1) =(¬NAD(t)) 
NAD(t + 1) =AMPK(t) 
Oxygen(t + 1) =Oxygen(t) 
SMAD6(t + 1) =CSL(t)
SMAD1(t + 1) =(¬SMAD6(t)) ∧ (¬NRP1(t)) ∧ ALK1(t) ∧ ENG(t) ∧ SMAD4(t)
SMAD2(t + 1) =(¬SMAD6(t)) ∧ (¬NRP1(t)) ∧ ALK5(t) ∧ BG(t) ∧ SMAD4(t)
Lactate(t + 1) =¬Oxygen(t)
PHDs(t + 1) =Oxygen(t) ∧ ¬Lactate(t)
HIF2a(t + 1) =SIRT1(t) ∨ (¬Oxygen(t))
LSD1(t + 1) =SIRT1(t) ∧ (¬NICD(t))
SIRT1(t + 1) =(HIF2a(t) ∨ HIF1(t) ∨ FOXO1(t)) ∧ NAD(t)
WNT5a(t + 1) =WNT5a(t)
ROR(t + 1) =WNT5a(t)
JNK(t + 1) =(ROR(t) ∨ RAC1(t)) ∧ DSH(t)
AP1(t + 1) =JNK(t)
FZD45(t + 1) =WNT5a(t)
WNT11(t + 1) =WNT11(t)
FZD6(t + 1) =WNT11(t) ∨ WNT5a(t)
RAC1(t + 1) =FZD6(t) ∧ DSH(t)
DAAM1(t + 1) =FZD6(t) ∨ DSH(t)
CDC42(t + 1) =DAAM1(t)
RHO(t + 1) =DAAM1(t) ∨ NRP1(t)
WNT7a(t + 1) =WNT7a(t) ∧ (¬sFRP1(t))
sFRP1(t + 1) =0
FZD47(t + 1) =WNT7a(t)
LRP56(t + 1) =WNT7a(t) ∧ (¬DKK1(t))
GSK3B(t + 1) =¬DSH(t) ∧ ¬LRP56(t)
Axin2(t + 1) =LEF1(t)
BTrCP(t + 1) =0 
Bcatenin(t + 1) =¬GSK3B(t) ∧ ¬BTrCP (t)
GROUCHO(t + 1) =¬Bcatenin(t) 
LEF1(t + 1) =(Bcatenin(t) ∧ ¬GROUCHO(t)) ∧ (NRARP(t) ∨ LEF1(t)) 
CyclinD1(t + 1) =Bcatenin(t) ∧ LEF1(t) 
FGF(t + 1) =FGF(t) 
FGFR2(t + 1) =FGF(t) 
FRS2a(t + 1) =FGFR2(t) 
JAGa(t + 1) =SMAD1(t) ∨ Bcatenin(t) 
JAGp(t + 1) =JAGp(t) 
DLL4a(t + 1) =ETS(t) ∨ CSL(t)
DKK1(t + 1) =0
DSH(t + 1) =(FZD45(t) ∧ ROR(t)) ∨ FZD6(t) ∨ (FZD47(t) ∧ LRP56(t))
DLL4p(t + 1) =DLL4p(t)
ADAM10(t + 1) =1
gSecretase(t + 1) =1
Notch4(t + 1) =ETS(t)
NOTCH(t + 1) =Notch4(t) ∧ gSecretase(t) ∧ ADAM10(t) ∧ DLL4p(t) ∧ (¬JAGp(t))
NICD(t + 1) =NOTCH(t) ∧ (¬NRARP(t))
CSL(t + 1) =NICD(t)
NRARP(t + 1) =CSL(t)
HEY1(t + 1) =CSL(t) ∨ ((SMAD1(t) ∨ SMAD2(t)) ∧ (¬SIRT1(t) ∧ ¬LSD1(t)))
Vegfr2(t + 1) =(ETS(t) ∧ ¬HEY1(t)) ∨ Lactate(t)
Vegfr1(t + 1) =CSL(t) ∨ ((HIF1(t) ∨ ETS(t)) ∧ (¬NFAT(t)))
PlGF(t + 1) =PlGF(t)
VEGFR1s(t + 1) =PlGF(t) ∧ Vegfr1(t)
VEGFB(t + 1) =VEGFB(t)
VEGFR11(t + 1) =Vegfr1(t)
VEGFR12(t + 1) =Vegfr1(t) ∧ Vegfr2(t)
IGF(t + 1) =IGF(t)
PKC(t + 1) =IGF(t)
ASFSF2(t+1)=PKC(t)
VEGFAxxxP(t + 1) =VEGFAxxxP(t)
VEGFAxxx(t + 1) =((¬(VEGFR1s(t) ∨ VEGFR12(t) ∨ VEGFR11(t)) ∧ VegfA(t) ∧ ART1(t) ∧ ASFSF2(t)) ∨ VEGFAxxxP(t))
p38MAPK(t + 1) =TAK1(t) ∨ SRC(t) ∨ DAG(t)
CLK1(t + 1) =p38MAPK(t)
CLK4(t + 1) =p38MAPK(t)
SRP55(t + 1) =CLK4(t) ∨ CLK1(t)
VEGFAxxxd(t + 1) =VegfA(t) ∧ SRP55(t)
VEGFCDp(t + 1) =VEGFCDp(t)
VEGFCD(t+1)=VEGFCD(t)
VEGFR22(t + 1) =Vegfr2(t) ∧ (PECAM1(t) ∨ ((VEGFCDp(t) ∨ VEGFAxxx(t)) ∧ ¬(VEGFAxxxd(t) ∨ VEPTP (t))))
VEGFR3(t + 1) =CSL(t)
VEGFR23(t + 1) =Vegfr2(t) ∧ VEGFR3(t) ∧ (PECAM1(t) ∨ VEGFAxxx(t) ∨ VEGFCDp(t))
VEGFR33(t + 1) =VEGFR3(t) ∧ (PECAM1(t) ∨ VEGFCDp(t) ∨ VEGFCD(t))
STAT3(t + 1) =VEGFR22(t)
VegfA(t + 1) =STAT3(t) ∨ NFAT(t) ∨ HIF1(t) ∨ Lactate(t) ∨ FOXO1(t) ∨ KLF2(t)
NRP1(t + 1) =(VEGFAxxx ∨ VEGFCDp) ∧ (¬CSL ∨ ETS)'''



a = a.replace(' ∨ ',' OR ').replace(' ∧ ',' AND ').replace('∧ ',' AND ').replace('¬',' NOT ').replace('(t)',' ').replace('(t + 1)',' ').replace('(t+1)',' ').replace('=',' = ').replace('  ',' ').replace('  ',' ').replace('  ',' ')

g = open('update_rules_models_in_literature_we_randomly_come_across/29230182.txt','w')
for line in a.split('\n'):
    if line.split(' = ')[1].strip() in [line.split(' = ')[0].strip(),'1','0']:
        print('detected constant: '+line)
    else:
        g.write(line+'\n')
g.close()
