const BaseRules = OrderedRuleSet(
    "Base",
    ("LowCounts", r -> get(r, :COUNTS, 2001.0) < 2000.0),
    (
        "Maraging",
        r -> (r.FIRSTELM == n"Fe") &&
                 (r.NI > 5) && (r.CO > 2) && sum(r, :FE, :NI, :CO, :MO) > 80.0,
    ),
    ("Quartz-like", r -> r.SI > 80),
    (
        "Biotite",
        r -> (r.K > 5) &&
                 (r.MG + r.FE > 20) &&
                 (r.AL > 10) &&
                 (r.SI > 10) && (r.K + r.MG + r.FE + r.AL + r.SI > 80) && (r.O > 5),
    ),
    (
        "Plagioclaise",
        r -> (r.CA + r.NA > 20) &&
                 (r.AL > 10) && (r.SI > 10) && (r.NA + r.CA + r.AL + r.SI > 80) && (r.O > 5),
    ),
    (
        "Othoclaise",
        r -> (r.K > 5) &&
                 (r.AL > 20) && (r.SI > 20) && (r.O > 10) && (r.K + r.AL + r.SI > 80),
    ),
    ("Dolomite", r -> (r.CA > 40) && (r.MG > 10) && (r.CA + r.MG > 70)),
    ("Olivine", r -> (r.MG + r.FE > 40) && (r.SI > 20) && (r.MG + r.FE + r.SI > 80)),
    (
        "Feldspar",
        r -> (r.K + r.NA + r.CA > 20) &&
                 (r.AL > 10) && (r.SI > 20) && (r.K + r.NA + r.CA + r.AL + r.SI > 80),
    ),
    ("Willemite", r -> (r.ZN > 20) && (r.SI > 20) && (r.ZN + r.SI > 80)),
    (
        "Garnet",
        r -> (r.MG + r.FE + r.MN + r.CA > 10) &&
                 (r.AL > 5) && (r.SI > 10) && (r.MG + r.FE + r.MN + r.CA + r.AL + r.SI > 80),
    ),
    (
        "Al-Cu alloy (2XXX)",
        r -> (r.AL > 60) && (r.SECONDELM == n"Cu") && (r.CU > 3) && (r.AL + r.CU > 90),
    ),
    (
        "Al-Mn alloy (3XXX)",
        r -> (r.AL > 60) && (r.SECONDELM == n"Mn") && (r.MN > 3) && (r.AL + r.MN > 90),
    ),
    (
        "Al-Si alloy (4XXX)",
        r -> (r.AL > 80) && (r.SECONDELM == n"Si") && (r.SI > 1) && (r.AL + r.SI > 90),
    ),
    (
        "Al-Mg alloy (5XXX)",
        r -> (r.FIRSTELM == n"Al") && (r.SECONDELM == n"Mg") && (r.MG > 5) && (r.AL + r.MG > 80),
    ),
    (
        "Al-Si-Mg (6XXX)",
        r -> (r.AL > 90) &&
                 (r.SI > 0) && (r.MG > 0) && (r.AL + r.SI + r.MG + r.MN + r.CU + r.FE > 98),
    ),
    (
        "Al-Zn alloy (7XXX)",
        r -> (r.AL > 80) && (r.SECONDELM == n"Zn") && (r.ZN > 2) && (r.AL + r.ZN > 90),
    ),
    ("Al-Fe alloy", r -> (r.AL > 80) && (r.SECONDELM == n"Fe") && (r.FE > 2) && (r.AL + r.FE > 90)),
    ("Aluminosilicate", r -> (r.AL > 30) && (r.SI > 10) && (r.AL + r.SI > 80)),
    ("Zircon", r -> (r.ZR > 20) && (r.ZR + r.SI > 80)),
    (
        "Elgiloy",
        r -> (r.CO > 20) &&
                 (r.CR > 10) &&
                 (r.NI > 5) &&
                 (r.FE > 5) &&
                 ((r.MO > 2) || (r.S > 2)) &&
                 (r.CO + r.CR + r.NI + r.FE + r.MO + r.MN + r.S > 80),
    ),
    ("Aluminum", r -> r.AL > 80),
    ("Fe-Al alloy", r -> (r.AL > 20) && (r.FE > 20) && (r.AL + r.FE > 80)),
    ("Calcium silicate", r -> (r.CA > 40) && (r.SI > 10) && (r.CA + r.SI > 80)),
    ("Fe-Cr stainless", r -> (r.FE > 50) && (r.CR > 5) && (r.NI < 3) && (r.FE + r.CR > 90)),
    (
        "Fe-Cr-Ni Stainless",
        r -> (r.FIRSTELM == n"Fe") && (r.CR > 5) && (r.NI > 5) && (r.FE + r.CR + r.NI > 80),
    ),
    ("Fe-Ni stainless", r -> (r.FE > 50) && (r.NI > 10) && (r.FE + r.NI > 90)),
    ("Iron oxide", r -> (r.FE > 80) && (r.O > 10)),
    ("Iron", r -> r.FE > 80),
    ("Calcite", r -> (r.CA > 80) && (r.O > 5)),
    ("Zinc", r -> r.ZN > 80),
    ("Silicate", r -> r.FIRSTELM == n"Si"),
    (
        "Anglesite",
        r -> (r.FIRSTELM == n"Pb") && (r.SECONDELM == n"S") && (r.PB + r.S > 80) && (r.O > 5),
    ),
    ("Hashemite", r -> (r.BA > 20) && (r.CR > 20) && (r.BA + r.CR > 80) && (r.O > 5)),
    ("Barite", r -> (r.BA > 40) && (r.S > 10) && (r.BA + r.S > 80) && (r.O > 5)),
    ("Gypsum", r -> (r.CA > 20) && (r.S > 20) && (r.CA + r.S > 80) && (r.O > 5)),
    ("Lead", r -> (r.PB > 60) && (r.PB + r.S > 80)),
    ("Salt", r -> (r.K + r.NA > 20) && (r.CL > 20) && (r.K + r.NA + r.CL > 80)),
    ("Bismuth", r -> (r.BI > 80)),
    ("Brass", r -> (r.CU > 20) && (r.ZN > 20) && (r.CU + r.ZN > 80)),
    ("Bronze", r -> (r.CU > 20) && (r.SN > 20) && (r.CU + r.SN > 80)),
    ("Cupronickel", r -> (r.CU > 30) && (r.NI > 10) && (r.CU + r.NI > 80)),
    ("Sodium sulfide", r -> (r.NA > 20) && (r.S > 10) && (r.NA + r.S > 80)),
    ("Wurtzite (ZnS)", r -> (r.ZN > 20) && (r.S > 10) && (r.ZN + r.S > 80)),
    ("Pyrite (FeS2)", r -> (r.FE > 20) && (r.S > 20) && (r.FE + r.S > 80)),
    ("Silver sulfide", r -> (r.AG > 20) && (r.S > 10) && (r.AG + r.S > 80)),
    (
        "Sn62 solder",
        r -> (r.SN > 30) && (r.PB > 10) && (r.AG > 1) && (r.SN + r.PB + r.AG > 80),
    ),
    (
        "SACZ solder",
        r -> (r.SN + r.AG + r.CU + r.ZN > 80) &&
                 (r.SN > 10) && (r.CU > 10) && (r.AG > 10) && (r.ZN > 2),
    ),
    ("Lead solder", r -> (r.PB > 20) && (r.SN > 20) && (r.PB + r.SN > 80)),
    (
        "SAC solder",
        r -> (r.SN + r.AG + r.CU > 80) && (r.SN > 10) && (r.CU > 10) && (r.AG > 10),
    ),
    ("Titanium", r -> (r.TI > 80)),
    ("Calcium phosphate", r -> (r.CA > 40) && (r.P > 10) && (r.CA + r.P > 80)),
    ("Calcium fluoride", r -> (r.CA > 40) && (r.F > 10) && (r.CA + r.F > 80)),
    ("Iron-bearing", r -> (r.FE > 50)),
    ("Ca-bearing", r -> r.CA > 50),
    ("Other Silicate", r -> (r.SI > 20) && (r.AL + r.FE + r.CA + r.SI > 80)),
    ("Gold-Silver alloy", r -> (r.FIRSTELM == n"Au") && (r.SECONDELM == n"Ag") && (r.AU + r.AG > 80)),
    ("Gold - other", r -> r.AU > 80),
    ("Tin", r -> r.SN > 80),
    ("Celestine", r -> (r.SR > 40) && (r.S > 10) && (r.SR + r.S > 80) && (r.O > 5)),
    ("Barite-Celestine", r -> (r.BA + r.SR > 40) && (r.BA + r.SR + r.S > 80) && (r.O > 5)),
    ("Cu-rich", r -> r.CU > 80),
    ("Cr-bearing", r -> r.CR > 50),
    ("Silver sulfide", r -> (r.AG > 40) && (r.S > 10) && (r.AG + r.S > 80)),
    ("Silver", r -> r.AG > 80),
    ("MoS2", r -> r.MO > 40 && r.S > 30 && (r.MO + r.S > 90)),
    ("Monazite", r -> sum(r, :LA, :CE, :ND, :TH, :Y) > 30 && (r.P > 10) && sum(r, :LA, :CE, :ND, :TH, :Y, :P, :SI) > 80 ),
    ("Lanthanide", r -> sum(r, :LA, :CE, :ND, :Y) > 60),
    ("Silver chloride", r -> (r.AG > 40) && (r.CL > 10)),
    ("Magnesium oxide", r -> (r.MG > 80) && (r.O > 10)),
    ("Nickel", r -> r.NI > 80),
    ("Ca-Fe silicate", r -> (r.CA > 20) && (r.FE > 20) && (r.SI > 5)),
    ("Ca-Ti silicate", r -> (r.CA > 20) && (r.TI > 20) && (r.SI > 10)),
    ("Al-Zr-Cl", r -> (r.AL > 20) && (r.ZR > 20) && (r.CL > 5)),
    ("Chlorine-rich", r -> r.CL > 80),
    ("Au-Cu alloy", r -> (r.AU > 40) && (r.CU > 10) && (r.AU + r.CU + r.NI > 80)),
    ("Uranium-bearing", r -> r.U > 20),
    ("Thorium-bearing", r -> r.TH > 20),
    ("Platinum", r -> r.PT > 80)
)