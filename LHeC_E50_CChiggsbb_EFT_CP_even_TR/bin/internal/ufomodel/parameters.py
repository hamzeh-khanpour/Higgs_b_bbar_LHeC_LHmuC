# This file was automatically created by FeynRules 2.3.36
# Mathematica version: 10.0 for Linux x86 (64-bit) (December 4, 2014)
# Date: Fri 30 Sep 2022 19:13:49



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
TR = Parameter(name = 'TR',
               nature = 'external',
               type = 'real',
               value = 0.0,
               texname = 'T_R',
               lhablock = 'ALPPARS',
               lhacode = [ 1 ])

TI = Parameter(name = 'TI',
               nature = 'external',
               type = 'real',
               value = 0.0,
               texname = 'T_i',
               lhablock = 'ALPPARS',
               lhacode = [ 2 ])

cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.227736,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymdo = Parameter(name = 'ymdo',
                 nature = 'external',
                 type = 'real',
                 value = 0.00504,
                 texname = '\\text{ymdo}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 1 ])

ymup = Parameter(name = 'ymup',
                 nature = 'external',
                 type = 'real',
                 value = 0.00255,
                 texname = '\\text{ymup}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 2 ])

yms = Parameter(name = 'yms',
                nature = 'external',
                type = 'real',
                value = 0.101,
                texname = '\\text{yms}',
                lhablock = 'YUKAWA',
                lhacode = [ 3 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.27,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

yme = Parameter(name = 'yme',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{yme}',
                lhablock = 'YUKAWA',
                lhacode = [ 11 ])

ymm = Parameter(name = 'ymm',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{ymm}',
                lhablock = 'YUKAWA',
                lhacode = [ 13 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MMU = Parameter(name = 'MMU',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{MMU}',
                lhablock = 'MASS',
                lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.00255,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.27,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.00504,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.101,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

MBbarMB = Parameter(name = 'MBbarMB',
                    nature = 'internal',
                    type = 'real',
                    value = '4.27285',
                    texname = '\\text{MBbar}_{\\text{MB}}')

Nf = Parameter(name = 'Nf',
               nature = 'internal',
               type = 'real',
               value = '5',
               texname = 'N_f')

lambdaQCD5 = Parameter(name = 'lambdaQCD5',
                       nature = 'internal',
                       type = 'real',
                       value = '0.216',
                       texname = '\\lambda _{\\text{QCD}}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM1x1}')

CKM1x2 = Parameter(name = 'CKM1x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.sin(cabi)',
                   texname = '\\text{CKM1x2}')

CKM1x3 = Parameter(name = 'CKM1x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM1x3}')

CKM2x1 = Parameter(name = 'CKM2x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '-cmath.sin(cabi)',
                   texname = '\\text{CKM2x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM2x2}')

CKM2x3 = Parameter(name = 'CKM2x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM2x3}')

CKM3x1 = Parameter(name = 'CKM3x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM3x1}')

CKM3x2 = Parameter(name = 'CKM3x2',
                   nature = 'internal',
                   type = 'complex',
                   value = '0',
                   texname = '\\text{CKM3x2}')

CKM3x3 = Parameter(name = 'CKM3x3',
                   nature = 'internal',
                   type = 'complex',
                   value = '1',
                   texname = '\\text{CKM3x3}')

beta0 = Parameter(name = 'beta0',
                  nature = 'internal',
                  type = 'real',
                  value = '11 - (2*Nf)/3.',
                  texname = '\\beta _0')

beta1 = Parameter(name = 'beta1',
                  nature = 'internal',
                  type = 'real',
                  value = '51 - (19*Nf)/3.',
                  texname = '\\beta _1')

beta2 = Parameter(name = 'beta2',
                  nature = 'internal',
                  type = 'real',
                  value = '2857 - (5033*Nf)/9. + (325*Nf**2)/27.',
                  texname = '\\beta _2')

lmuMB = Parameter(name = 'lmuMB',
                  nature = 'internal',
                  type = 'real',
                  value = 'cmath.log(MB**2/lambdaQCD5**2)',
                  texname = 'l_{\\text{muMB}}')

lmuMH = Parameter(name = 'lmuMH',
                  nature = 'internal',
                  type = 'real',
                  value = 'cmath.log(MH**2/lambdaQCD5**2)',
                  texname = 'l_{\\text{muMH}}')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

alphaSbarMB = Parameter(name = 'alphaSbarMB',
                        nature = 'internal',
                        type = 'real',
                        value = '(4*cmath.pi*(1 + (4*beta1**2*(-1.25 + (beta0*beta2)/(8.*beta1**2) + (-0.5 + cmath.log(lmuMB))**2))/(beta0**4*lmuMB**2) - (2*beta1*cmath.log(lmuMB))/(beta0**2*lmuMB)))/(beta0*lmuMB)',
                        texname = '\\alpha _{\\text{sMB}}')

alphaSbarMH = Parameter(name = 'alphaSbarMH',
                        nature = 'internal',
                        type = 'real',
                        value = '(4*cmath.pi*(1 + (4*beta1**2*(-1.25 + (beta0*beta2)/(8.*beta1**2) + (-0.5 + cmath.log(lmuMH))**2))/(beta0**4*lmuMH**2) - (2*beta1*cmath.log(lmuMH))/(beta0**2*lmuMH)))/(beta0*lmuMH)',
                        texname = '\\alpha _{\\text{sMH}}')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

xc1b = Parameter(name = 'xc1b',
                 nature = 'internal',
                 type = 'real',
                 value = 'alphaSbarMH/cmath.pi',
                 texname = 'x_{\\text{c1}_b}')

xc2b = Parameter(name = 'xc2b',
                 nature = 'internal',
                 type = 'real',
                 value = 'alphaSbarMB/cmath.pi',
                 texname = 'x_{\\text{c2}_b}')

c1b = Parameter(name = 'c1b',
                nature = 'internal',
                type = 'real',
                value = '2.0159290713817337*xc1b**0.52174*(1 + 1.175*xc1b + 1.501*xc1b**2 + 0.1725*xc1b**3)',
                texname = '\\text{c1}_b')

c2b = Parameter(name = 'c2b',
                nature = 'internal',
                type = 'real',
                value = '2.0159290713817337*xc2b**0.52174*(1 + 1.175*xc2b + 1.501*xc2b**2 + 0.1725*xc2b**3)',
                texname = '\\text{c2}_b')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

MBbarMH = Parameter(name = 'MBbarMH',
                    nature = 'internal',
                    type = 'real',
                    value = '(c1b*MBbarMB)/c2b',
                    texname = '\\text{MBbar}_{\\text{MH}}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2))/vev',
               texname = '\\text{yc}')

ydo = Parameter(name = 'ydo',
                nature = 'internal',
                type = 'real',
                value = '(ymdo*cmath.sqrt(2))/vev',
                texname = '\\text{ydo}')

ye = Parameter(name = 'ye',
               nature = 'internal',
               type = 'real',
               value = '(yme*cmath.sqrt(2))/vev',
               texname = '\\text{ye}')

ym = Parameter(name = 'ym',
               nature = 'internal',
               type = 'real',
               value = '(ymm*cmath.sqrt(2))/vev',
               texname = '\\text{ym}')

ys = Parameter(name = 'ys',
               nature = 'internal',
               type = 'real',
               value = '(yms*cmath.sqrt(2))/vev',
               texname = '\\text{ys}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

yup = Parameter(name = 'yup',
                nature = 'internal',
                type = 'real',
                value = '(ymup*cmath.sqrt(2))/vev',
                texname = '\\text{yup}')

deltaQCDb = Parameter(name = 'deltaQCDb',
                      nature = 'internal',
                      type = 'real',
                      value = '(17*alphaSbarMH)/(3.*cmath.pi) + (alphaSbarMH**2*(35.94 - 1.36*Nf))/cmath.pi**2 + (alphaSbarMH**2*(1.57 + cmath.log(MBbarMH**2/MH**2)/9. - (2*cmath.log(MH**2/MT**2))/3.))/cmath.pi**2',
                      texname = '\\delta _{\\text{QCDb}}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

widthHbb = Parameter(name = 'widthHbb',
                     nature = 'internal',
                     type = 'real',
                     value = '(3*(1 + deltaQCDb)*Gf*MBbarMH**2*MH*cmath.sqrt(1 - (4*MB**2)/MH**2))/(4.*cmath.pi*cmath.sqrt(2))',
                     texname = '\\text{width}_{\\text{bb}}')

widthHbbBSM = Parameter(name = 'widthHbbBSM',
                        nature = 'internal',
                        type = 'real',
                        value = '(9*TI**2 + 6*TR + 9*TR**2)*widthHbb',
                        texname = '\\text{width}_{\\text{bb}}')

WH = Parameter(name = 'WH',
               nature = 'internal',
               type = 'real',
               value = '0.00412 + widthHbbBSM',
               texname = 'W_{\\text{Higgs}}')

