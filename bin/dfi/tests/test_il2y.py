import pandas as pd
import numpy as np
import dfi
from dfi.datafiles import example_covar, example_pdb
from dfi.dfi_calc import check_args
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def test_Il2Y():
    sample = StringIO("""
    ResI,ChainID,Res,R,dfi,pctdfi,fdfi,pctfdfi,adfi,ravg,rmin,A
1,A,ASN,N,0.0739300399161,0.95,0.733721,0.55,0.45,14.8066793711,14.8066793711,NotA
2,A,LEU,L,0.0570598728478,0.75,0.760436,0.65,0.5,13.3382425379,13.3382425379,NotA
3,A,TYR,Y,0.0366594203159,0.3,0.767825,0.7,0.45,10.9987171979,10.9987171979,NotA
4,A,ILE,I,0.0264925229724,0.15,0.887541,0.85,0.45,10.0844360279,10.0844360279,A
5,A,GLN,Q,0.0313964597566,0.25,0.71392,0.4,-0.05,9.31297057871,9.31297057871,NotA
6,A,TRP,W,0.0219462741489,0.05,0.656116,0.2,-0.4,7.08785870344,7.08785870344,NotA
7,A,LEU,L,0.0301269473082,0.2,0.73717,0.6,-0.15,5.2364966342,5.2364966342,NotA
8,A,LYS,K,0.0470417632665,0.65,0.697586,0.25,-0.45,5.47704564524,5.47704564524,NotA
9,A,ASP,D,0.0460315285765,0.6,0.959495,0.9,0.0,3.89091814357,3.89091814357,NotA
10,A,GLY,G,0.0379863773561,0.4,1.528388,1.0,0.0,0.0,0.0,NotA
11,A,GLY,G,0.0259132562681,0.1,1.182518,0.95,0.0,3.8635549692,3.8635549692,NotA
12,A,PRO,P,0.0398701906642,0.55,0.712332,0.35,-0.3,6.31214092048,6.31214092048,NotA
13,A,SER,S,0.0383509211918,0.5,0.413331,0.05,-0.8,5.13002563346,5.13002563346,NotA
14,A,SER,S,0.0380024946005,0.45,0.551883,0.1,-0.7,5.21916794135,5.21916794135,NotA
15,A,GLY,G,0.0548582188702,0.7,0.824349,0.75,0.2,8.88162451357,8.88162451357,NotA
16,A,ARG,R,0.0609643869267,0.85,0.718268,0.45,-0.05,9.24985697187,9.24985697187,NotA
17,A,PRO,P,0.0727215382269,0.9,0.725202,0.5,0.2,10.8655981888,10.8655981888,NotA
18,A,PRO,P,0.0372565419539,0.35,0.570706,0.15,-0.2,10.7198590942,10.7198590942,NotA
19,A,PRO,P,0.0588722376646,0.8,0.70863,0.3,0.1,13.3328339823,13.3328339823,NotA
20,A,SER,S,0.164519007168,1.0,0.878281,0.8,0.75,16.6916192444,16.6916192444,A
""")

    df = pd.read_csv(sample)
    print(df)
    sysls = ['--pdb', example_pdb, '--fdfi', 'A10']
    pdbfile, pdbid, mdhess, ls_reschain, chain_nam = check_args(
        sysls)
    df_dfi = dfi.dfi_calc.calc_dfi(
        pdbfile, covar=mdhess, ls_reschain=ls_reschain, chain_name=chain_nam)
    assert np.all(df_dfi.Res.values == df.Res.values)
#    assert np.allclose(df_dfi.pctdfi.values, df.pctdfi.values)
#    assert np.allclose(df_dfi.pctfdfi.values, df.pctfdfi.values)
#    assert np.allclose(df_dfi.ravg.values, df.ravg.values)

    sysls = ['--pdb', example_pdb, '--covar', example_covar,
             '--fdfi', 'A10']
    pdbfile, pdbid, mdhess, ls_reschain, chain_nam = check_args(
        sysls)
    df_dfi = dfi.dfi_calc.calc_dfi(
        pdbfile, covar=mdhess, ls_reschain=ls_reschain, chain_name=chain_nam)
#    assert np.all(df_dfi.Res.values == df.Res.values)
#    assert np.allclose(df_dfi.pctdfi.values, df.pctdfi.values)
#    assert np.allclose(df_dfi.pctfdfi.values, df.pctfdfi.values)
#    assert np.allclose(df_dfi.ravg.values, df.ravg.values)
