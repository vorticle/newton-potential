/*
 * Copyright (C) 2020 Matthias Kirchhart
 *
 * This is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this software; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#include <geometry/quadrules.h>

namespace
{

using geometry::interval_quad_node;

const std::array<interval_quad_node,1> legendre_rule_01
{{
    interval_quad_node { 0.5, 1 }
}};

const std::array<interval_quad_node,2> legendre_rule_02
{{
    interval_quad_node { 0.21132486540518711774542560974902, 0.5 },
    interval_quad_node { 0.78867513459481288225457439025098, 0.5 }
}};

const std::array<interval_quad_node,3> legendre_rule_03
{{
    interval_quad_node { 0.5,0.44444444444444444444444444444444 },
    interval_quad_node { 0.11270166537925831148207346002176, 0.27777777777777777777777777777778 },
    interval_quad_node { 0.88729833462074168851792653997824, 0.27777777777777777777777777777778 }
}};

const std::array<interval_quad_node,4> legendre_rule_04
{{
    interval_quad_node { 0.33000947820757186759866712044838,  0.326072577431273071313468025389 },
    interval_quad_node { 0.66999052179242813240133287955162,  0.326072577431273071313468025389 },
    interval_quad_node { 0.069431844202973712388026755553595, 0.173927422568726928686531974611 },
    interval_quad_node { 0.9305681557970262876119732444464,   0.173927422568726928686531974611 }
}};

const std::array<interval_quad_node,5> legendre_rule_05
{{
    interval_quad_node { 0.5,0.28444444444444444444444444444444 },
    interval_quad_node { 0.2307653449471584544818427896499,   0.23931433524968323402064575741782 },
    interval_quad_node { 0.7692346550528415455181572103501,   0.23931433524968323402064575741782 },
    interval_quad_node { 0.046910077030668003601186560850304, 0.11846344252809454375713202035996 },
    interval_quad_node { 0.9530899229693319963988134391497,   0.11846344252809454375713202035996 }
}};

const std::array<interval_quad_node,6> legendre_rule_06
{{
    interval_quad_node { 0.83060469323313225683069979750995,  0.18038078652406930378491675691886  },
    interval_quad_node { 0.16939530676686774316930020249005,  0.18038078652406930378491675691886  },
    interval_quad_node { 0.38069040695840154568474913915964,  0.23395696728634552369493517199478  },
    interval_quad_node { 0.61930959304159845431525086084036,  0.23395696728634552369493517199478  },
    interval_quad_node { 0.033765242898423986093849222753003, 0.085662246189585172520148071086366 },
    interval_quad_node { 0.966234757101576013906150777247,    0.085662246189585172520148071086366 }
}};

const std::array<interval_quad_node,7> legendre_rule_07
{{
    interval_quad_node { 0.5, 0.2089795918367346938775510204081 },
    interval_quad_node { 0.70292257568869858345330320603848,  0.19091502525255947247518488774449  },
    interval_quad_node { 0.29707742431130141654669679396152,  0.19091502525255947247518488774449  },
    interval_quad_node { 0.12923440720030278006806761335961,  0.13985269574463833395073388571189  },
    interval_quad_node { 0.87076559279969721993193238664039,  0.13985269574463833395073388571189  },
    interval_quad_node { 0.025446043828620737736905157976074, 0.064742483084434846635305716339541 },
    interval_quad_node { 0.97455395617137926226309484202393,  0.064742483084434846635305716339541 }
}};

const std::array<interval_quad_node,8> legendre_rule_08
{{
    interval_quad_node { 0.40828267875217509753026192881991,  0.1813418916891809914825752246386   },
    interval_quad_node { 0.59171732124782490246973807118009,  0.1813418916891809914825752246386   },
    interval_quad_node { 0.23723379504183550709113047540538,  0.1568533229389436436689811009933   },
    interval_quad_node { 0.76276620495816449290886952459462,  0.1568533229389436436689811009933   },
    interval_quad_node { 0.10166676129318663020422303176208,  0.11119051722668723527217799721312  },
    interval_quad_node { 0.89833323870681336979577696823792,  0.11119051722668723527217799721312  },
    interval_quad_node { 0.019855071751231884158219565715264, 0.050614268145188129576265677154981 },
    interval_quad_node { 0.98014492824876811584178043428474,  0.050614268145188129576265677154981 }
}};

const std::array<interval_quad_node,9> legendre_rule_09
{{
    interval_quad_node { 0.5,0.16511967750062988158226253464349 },
    interval_quad_node { 0.081984446336682102850285105965133,0.090324080347428702029236015621456 },
    interval_quad_node { 0.91801555366331789714971489403487,0.090324080347428702029236015621456  },
    interval_quad_node { 0.015919880246186955082211898548164,0.040637194180787205985946079055262 },
    interval_quad_node { 0.98408011975381304491778810145184,0.040637194180787205985946079055262  },
    interval_quad_node { 0.33787328829809553548073099267833,0.15617353852000142003431520329222   },
    interval_quad_node { 0.66212671170190446451926900732167,0.15617353852000142003431520329222   },
    interval_quad_node { 0.19331428364970480134564898032926,0.13030534820146773115937143470932   },
    interval_quad_node { 0.80668571635029519865435101967074,0.13030534820146773115937143470932   }
}};

const std::array<interval_quad_node,10> legendre_rule_10
{{
    interval_quad_node { 0.42556283050918439455758699943514,0.14776211235737643508694649732567   },
    interval_quad_node { 0.57443716949081560544241300056486,0.14776211235737643508694649732567   },
    interval_quad_node { 0.28330230293537640460036702841711,0.13463335965499817754561346078473   },
    interval_quad_node { 0.71669769706462359539963297158289,0.13463335965499817754561346078473   },
    interval_quad_node { 0.16029521585048779688283631744256,0.10954318125799102199776746711408   },
    interval_quad_node { 0.83970478414951220311716368255744,0.10954318125799102199776746711408   },
    interval_quad_node { 0.067468316655507744633951655788253,0.074725674575290296572888169828849 },
    interval_quad_node { 0.93253168334449225536604834421175,0.074725674575290296572888169828849  },
    interval_quad_node { 0.013046735741414139961017993957774,0.033335672154344068796784404946666 },
    interval_quad_node { 0.98695326425858586003898200604223,0.033335672154344068796784404946666  }
}};

const std::array<interval_quad_node,11> legendre_rule_11
{{
    interval_quad_node {0.50000000000000000000000000000000,0.13646254338895031535724176416817},
    interval_quad_node {0.36522842202382751383423400729957,0.13140227225512333109034443494525},
    interval_quad_node {0.63477157797617248616576599270043,0.13140227225512333109034443494525},
    interval_quad_node {0.24045193539659409203713716527070,0.11659688229599523995926185242159},
    interval_quad_node {0.75954806460340590796286283472930,0.11659688229599523995926185242159},
    interval_quad_node {0.13492399721297533795329187398442,0.093145105463867125713048820715828},
    interval_quad_node {0.86507600278702466204670812601558,0.093145105463867125713048820715828},
    interval_quad_node {0.056468700115952350462421115348036,0.062790184732452312317347149611970},
    interval_quad_node {0.94353129988404764953757888465196,0.062790184732452312317347149611970},
    interval_quad_node {0.010885670926971503598030999438571,0.027834283558086833241376860221274},
    interval_quad_node {0.98911432907302849640196900056143,0.027834283558086833241376860221274}
}};

const std::array<interval_quad_node,12> legendre_rule_12
{{
    interval_quad_node {0.43738329574426554226377931526807,0.12457352290670139250028121802148},
    interval_quad_node {0.56261670425573445773622068473193,0.12457352290670139250028121802148},
    interval_quad_node {0.31608425050090990312365423167814,0.11674626826917740438042494946244},
    interval_quad_node {0.68391574949909009687634576832186,0.11674626826917740438042494946244},
    interval_quad_node {0.20634102285669127635164879052973,0.10158371336153296087453222790490},
    interval_quad_node {0.79365897714330872364835120947027,0.10158371336153296087453222790490},
    interval_quad_node {0.11504866290284765648155308339359,0.080039164271673113167326264771680},
    interval_quad_node {0.88495133709715234351844691660641,0.080039164271673113167326264771680},
    interval_quad_node {0.047941371814762571660767066940452,0.053469662997659215480127359096998},
    interval_quad_node {0.95205862818523742833923293305955,0.053469662997659215480127359096998},
    interval_quad_node {0.0092196828766403746547254549253596,0.023587668193255913597307980742509},
    interval_quad_node {0.99078031712335962534527454507464,0.023587668193255913597307980742509}    
}};

const std::array<interval_quad_node,13> legendre_rule_13
{{
    interval_quad_node {0.50000000000000000000000000000000,0.11627577661543695509729475763442},
    interval_quad_node {0.38477084202243260296723593945101,0.11314159013144861920604509301989},
    interval_quad_node {0.61522915797756739703276406054899,0.11314159013144861920604509301989},
    interval_quad_node {0.27575362448177657356104357393618,0.10390802376844425115626160965303},
    interval_quad_node {0.72424637551822342643895642606382,0.10390802376844425115626160965303},
    interval_quad_node {0.17882533027982988967800769650224,0.089072990380972869140023345998049},
    interval_quad_node {0.82117466972017011032199230349776,0.089072990380972869140023345998049},
    interval_quad_node {0.099210954633345043602896755208570,0.069436755109893619231800888434436},
    interval_quad_node {0.90078904536665495639710324479143,0.069436755109893619231800888434436},
    interval_quad_node {0.041200800388511017396726081749640,0.046060749918864223957210887976899},
    interval_quad_node {0.95879919961148898260327391825036,0.046060749918864223957210887976899},
    interval_quad_node {0.0079084726407059252635852755964452,0.020242002382657939760010796100493},
    interval_quad_node {0.99209152735929407473641472440355,0.020242002382657939760010796100493}
}};

const std::array<interval_quad_node,14> legendre_rule_14
{{
    interval_quad_node {0.44597252564632816896687767489008,0.10763192673157889509793822165813},
    interval_quad_node {0.55402747435367183103312232510992,0.10763192673157889509793822165813},
    interval_quad_node {0.34044381553605511978216408791576,0.10259923186064780198296203283061},
    interval_quad_node {0.65955618446394488021783591208424,0.10259923186064780198296203283061},
    interval_quad_node {0.24237568182092295401735464072441,0.092769198738968906870858295062579},
    interval_quad_node {0.75762431817907704598264535927559,0.092769198738968906870858295062579},
    interval_quad_node {0.15635354759415726492599009849033,0.078601583579096767284800969311921},
    interval_quad_node {0.84364645240584273507400990150967,0.078601583579096767284800969311921},
    interval_quad_node {0.086399342465117503405102628674803,0.060759285343951592344707404536238},
    interval_quad_node {0.91360065753488249659489737132520,0.060759285343951592344707404536238},
    interval_quad_node {0.035782558168213241331804430311063,0.040079043579880104902816638531427},
    interval_quad_node {0.96421744183178675866819556968894,0.040079043579880104902816638531427},
    interval_quad_node {0.0068580956515938305792013666479736,0.017559730165875931515916438069096},
    interval_quad_node {0.99314190434840616942079863335203,0.017559730165875931515916438069096}
}};

const std::array<interval_quad_node,15> legendre_rule_15
{{
    interval_quad_node {0.50000000000000000000000000000000,0.10128912096278063644031009998376},
    interval_quad_node {0.39940295300128273884968584830270,0.099215742663555788228059163221920},
    interval_quad_node {0.60059704699871726115031415169730,0.099215742663555788228059163221920},
    interval_quad_node {0.30292432646121831505139631450948,0.093080500007781105513400280933211},
    interval_quad_node {0.69707567353878168494860368549052,0.093080500007781105513400280933211},
    interval_quad_node {0.21451391369573057623138663137304,0.083134602908496966776600430240604},
    interval_quad_node {0.78548608630426942376861336862696,0.083134602908496966776600430240604},
    interval_quad_node {0.13779113431991497629190697269303,0.069785338963077157223902397255514},
    interval_quad_node {0.86220886568008502370809302730697,0.069785338963077157223902397255514},
    interval_quad_node {0.075896708294786391899675839612892,0.053579610233585967505934773342935},
    interval_quad_node {0.92410329170521360810032416038711,0.053579610233585967505934773342935},
    interval_quad_node {0.031363303799647047846120526144895,0.035183023744054062354633708225334},
    interval_quad_node {0.96863669620035295215387947385510,0.035183023744054062354633708225334},
    interval_quad_node {0.0060037409897572857552171407066937,0.015376620998058634177314196788602},
    interval_quad_node {0.99399625901024271424478285929331,0.015376620998058634177314196788602}
}};

const std::array<interval_quad_node,16> legendre_rule_16
{{
    interval_quad_node {0.45249374508118127990734033228752,0.094725305227534248142698361604142},
    interval_quad_node {0.54750625491881872009265966771248,0.094725305227534248142698361604142},
    interval_quad_node {0.35919822461037054338476974926975,0.091301707522461794433381833984610},
    interval_quad_node {0.64080177538962945661523025073025,0.091301707522461794433381833984610},
    interval_quad_node {0.27099161117138630682879027850821,0.084578259697501269094656039515180},
    interval_quad_node {0.72900838882861369317120972149179,0.084578259697501269094656039515180},
    interval_quad_node {0.19106187779867812577666411797560,0.074797994408288366040750865273739},
    interval_quad_node {0.80893812220132187422333588202440,0.074797994408288366040750865273739},
    interval_quad_node {0.12229779582249848305244940257628,0.062314485627766936026238141096008},
    interval_quad_node {0.87770220417750151694755059742372,0.062314485627766936026238141096008},   
    interval_quad_node {0.067184398806084128059766051143803,0.047579255841246392404962553801123},
    interval_quad_node {0.93281560119391587194023394885620,0.047579255841246392404962553801123},
    interval_quad_node {0.027712488463383711961005792232696,0.031126761969323946431421918497189},
    interval_quad_node {0.97228751153661628803899420776730,0.031126761969323946431421918497189},
    interval_quad_node {0.0052995325041750337019229132748337,0.013576229705877047425890286228009},
    interval_quad_node {0.99470046749582496629807708672517,0.013576229705877047425890286228009}    
}};

const std::array<interval_quad_node,17> legendre_rule_17
{{
    interval_quad_node {0.50000000000000000000000000000000,0.089723235178103262729132822130943},
    interval_quad_node {0.41075790925207607207466125317297,0.088281352683496323162635495056599},
    interval_quad_node {0.58924209074792392792533874682703,0.088281352683496323162635495056599},
    interval_quad_node {0.32438411827306184235140724145233,0.084002051078225022254985331894162},
    interval_quad_node {0.67561588172693815764859275854767,0.084002051078225022254985331894162},
    interval_quad_node {0.24365473145676151605687671568522,0.077022880538405144040715797400979},
    interval_quad_node {0.75634526854323848394312328431478,0.077022880538405144040715797400979},
    interval_quad_node {0.17116442039165461707484889167850,0.067568184234262736643159990851175},
    interval_quad_node {0.82883557960834538292515110832150,0.067568184234262736643159990851175},
    interval_quad_node {0.10924299805159929653738497223976,0.055941923596701985547394192813178},
    interval_quad_node {0.89075700194840070346261502776024,0.055941923596701985547394192813178},
    interval_quad_node {0.059880423136507048938522152755922,0.042518074158589590441767685095531},
    interval_quad_node {0.94011957686349295106147784724408,0.042518074158589590441767685095531},
    interval_quad_node {0.024662239115616119388641521052098,0.027729764686993600564720082679122},
    interval_quad_node {0.97533776088438388061135847894790,0.027729764686993600564720082679122},
    interval_quad_node {0.0047122623427913321622829900296674,0.012074151434273965980055013143783},
    interval_quad_node {0.99528773765720866783771700997033,0.012074151434273965980055013143783}
}};

const std::array<interval_quad_node,18> legendre_rule_18
{{
    interval_quad_node {0.45761249347913234937886907353211,0.084571191481571795920328235067493},
    interval_quad_node {0.54238750652086765062113092646789,0.084571191481571795920328235067493},
    interval_quad_node {0.37405688715424724520551357256104,0.082138241872916361493026888232964},
    interval_quad_node {0.62594311284575275479448642743896,0.082138241872916361493026888232964},
    interval_quad_node {0.29412441926857867698203410308347,0.077342337563132622462709001918187},
    interval_quad_node {0.70587558073142132301796589691653,0.077342337563132622462709001918187},
    interval_quad_node {0.22011458446302623269606422573734,0.070321457335325325602365651875974},
    interval_quad_node {0.77988541553697376730393577426266,0.070321457335325325602365651875974},
    interval_quad_node {0.15415647846982339606255445935558,0.061277603355739230092259563400101},
    interval_quad_node {0.84584352153017660393744554064442,0.061277603355739230092259563400101},
    interval_quad_node {0.098147520513738442158791272492705,0.050471022053143582781406992462417},
    interval_quad_node {0.90185247948626155784120872750730,0.050471022053143582781406992462417},
    interval_quad_node {0.053698766751222130396969704436427,0.038212865127444528264564838808318},
    interval_quad_node {0.94630123324877786960303029556357,0.038212865127444528264564838808318},
    interval_quad_node {0.022088025214301122409402053535112,0.024857274447484898226667473101319},
    interval_quad_node {0.97791197478569887759059794646489,0.024857274447484898226667473101319},
    interval_quad_node {0.0042174157895345266349919976469246,0.010808006763241655156671355133226},
    interval_quad_node {0.99578258421046547336500800235308,0.010808006763241655156671355133226}
}};

const std::array<interval_quad_node,19> legendre_rule_19
{{
    interval_quad_node {0.50000000000000000000000000000000,0.080527224924391847989581812660458},
    interval_quad_node {0.41982067717988731206595194212963,0.079484421696977173824978219732524},
    interval_quad_node {0.58017932282011268793404805787037,0.079484421696977173824978219732524},
    interval_quad_node {0.34171795001818508400494133557508,0.076383021032929833389427700448831},
    interval_quad_node {0.65828204998181491599505866442492,0.076383021032929833389427700448831},
    interval_quad_node {0.26771462931201952714136642594795,0.071303351086803305887873054720951},
    interval_quad_node {0.73228537068798047285863357405205,0.071303351086803305887873054720951},
    interval_quad_node {0.19972734766915948826518091752688,0.064376981269668113837757892428439},
    interval_quad_node {0.80027265233084051173481908247312,0.064376981269668113837757892428439},
    interval_quad_node {0.13951691133238531069145206958811,0.055783322773666997358011950840883},
    interval_quad_node {0.86048308866761468930854793041189,0.055783322773666997358011950840883},
    interval_quad_node {0.088642671731428587510538756643643,0.045745010811224999732231047061920},
    interval_quad_node {0.91135732826857141248946124335636,0.045745010811224999732231047061920},
    interval_quad_node {0.048422048192591049178669535733844,0.034522271368820613290354129003007},
    interval_quad_node {0.95157795180740895082133046426616,0.034522271368820613290354129003007},
    interval_quad_node {0.019895923932584984573610579656174,0.022407113382849800166419078700997},
    interval_quad_node {0.98010407606741501542638942034383,0.022407113382849800166419078700997},
    interval_quad_node {0.0037965780782077984054911648733698,0.0097308941148632385181560207322192},
    interval_quad_node {0.99620342192179220159450883512663,0.0097308941148632385181560207322192}
}};

const std::array<interval_quad_node,20> legendre_rule_20
{{
    interval_quad_node {0.46173673943325133312267979530058,0.076376693565362925349042165977549},
    interval_quad_node {0.53826326056674866687732020469942,0.076376693565362925349042165977549},
    interval_quad_node {0.38610707442917746095975190231571,0.074586493236301873393914368500985},
    interval_quad_node {0.61389292557082253904024809768429,0.074586493236301873393914368500985},
    interval_quad_node {0.31314695564229021966372591148754,0.071048054659191025664649162533582},
    interval_quad_node {0.68685304435770978033627408851246,0.071048054659191025664649162533582},
    interval_quad_node {0.24456649902458645099781797452237,0.065844319224588313449247249874082},
    interval_quad_node {0.75543350097541354900218202547763,0.065844319224588313449247249874082},
    interval_quad_node {0.18197315963674248727358165188686,0.059097265980759208656188688855691},
    interval_quad_node {0.81802684036325751272641834811314,0.059097265980759208656188688855691},
    interval_quad_node {0.12683404676992460369284746482218,0.050965059908620217518375067740175},
    interval_quad_node {0.87316595323007539630715253517782,0.050965059908620217518375067740175},
    interval_quad_node {0.080441514088890588302735469149240,0.041638370788352374362379071611023},
    interval_quad_node {0.91955848591110941169726453085076,0.041638370788352374362379071611023},
    interval_quad_node {0.043882785874337047066123779398351,0.031336024167054531784753267593521},
    interval_quad_node {0.95611721412566295293387622060165,0.031336024167054531784753267593521},
    interval_quad_node {0.018014036361043104366166934401361,0.020300714900193470665519976137466},
    interval_quad_node {0.98198596363895689563383306559864,0.020300714900193470665519976137466},
    interval_quad_node {0.0034357004074525376069388057643399,0.0088070035695760591559309811759264},
    interval_quad_node {0.99656429959254746239306119423566,0.0088070035695760591559309811759264}
}};

} 

namespace geometry
{

const interval_quadrule get_interval_quadrule( size_t degree )
{
    size_t N;
    const interval_quad_node *rule;

    switch ( degree )
    {
    case  0: rule = legendre_rule_01.data(); N =  1; break;
    case  1: rule = legendre_rule_01.data(); N =  1; break;
    case  2: rule = legendre_rule_02.data(); N =  2; break;
    case  3: rule = legendre_rule_02.data(); N =  2; break;
    case  4: rule = legendre_rule_03.data(); N =  3; break;
    case  5: rule = legendre_rule_03.data(); N =  3; break;
    case  6: rule = legendre_rule_04.data(); N =  4; break;
    case  7: rule = legendre_rule_04.data(); N =  4; break;
    case  8: rule = legendre_rule_05.data(); N =  5; break;
    case  9: rule = legendre_rule_05.data(); N =  5; break;
    case 10: rule = legendre_rule_06.data(); N =  6; break;
    case 11: rule = legendre_rule_06.data(); N =  6; break;
    case 12: rule = legendre_rule_07.data(); N =  7; break;
    case 13: rule = legendre_rule_07.data(); N =  7; break;
    case 14: rule = legendre_rule_08.data(); N =  8; break;
    case 15: rule = legendre_rule_08.data(); N =  8; break;
    case 16: rule = legendre_rule_09.data(); N =  9; break;
    case 17: rule = legendre_rule_09.data(); N =  9; break;
    case 18: rule = legendre_rule_10.data(); N = 10; break;
    case 19: rule = legendre_rule_10.data(); N = 10; break;
    case 20: rule = legendre_rule_11.data(); N = 11; break;
    case 21: rule = legendre_rule_11.data(); N = 11; break;
    case 22: rule = legendre_rule_12.data(); N = 12; break;
    case 23: rule = legendre_rule_12.data(); N = 12; break;
    case 24: rule = legendre_rule_13.data(); N = 13; break;
    case 25: rule = legendre_rule_13.data(); N = 13; break;
    case 26: rule = legendre_rule_14.data(); N = 14; break;
    case 27: rule = legendre_rule_14.data(); N = 14; break;
    case 28: rule = legendre_rule_15.data(); N = 15; break;
    case 29: rule = legendre_rule_15.data(); N = 15; break;
    case 30: rule = legendre_rule_16.data(); N = 16; break;
    case 31: rule = legendre_rule_16.data(); N = 16; break;
    case 32: rule = legendre_rule_17.data(); N = 17; break;
    case 33: rule = legendre_rule_17.data(); N = 17; break;
    case 34: rule = legendre_rule_18.data(); N = 18; break;
    case 35: rule = legendre_rule_18.data(); N = 18; break;
    case 36: rule = legendre_rule_19.data(); N = 19; break;
    case 37: rule = legendre_rule_19.data(); N = 19; break;
    case 38: rule = legendre_rule_20.data(); N = 20; break;
    case 39: rule = legendre_rule_20.data(); N = 20; break;
    default:
        throw std::logic_error { "Requested too high order interval quadrature rule." }; 
    }

    interval_quadrule result( rule, rule + N );
    return result;
}

interval_quadrule get_refined_interval_quadrule( size_t degree, size_t level )
{
    const interval_quadrule rule( get_interval_quadrule(degree) );
    size_t N { rule.size() };

    size_t M   { 1ul << level };
    real   dx  { 1./M };
    interval_quadrule result( N*M );
    for ( size_t i = 0; i < M; ++i )
    for ( size_t j = 0; j < N; ++j )
    {
        result[ N*i + j ].x = dx*(i+rule[ j ].x);
        result[ N*i + j ].w = dx*rule[ j ].w;
    }

    return result;
}

}
