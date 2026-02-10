{\rtf1\ansi\ansicpg936\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # \uc0\u26680 \u24515 \u26399 \u21002 \u21015 \u34920  (JMCB - CNS \u32423 \u21035 )\
# \uc0\u21482 \u35201 \u26399 \u21002 \u21517 \u21253 \u21547 \u20197 \u19979 \u23383 \u31526 \u20018 \u65292 \u23601 \u20250 \u34987 \u20445 \u30041 \
TARGET_JOURNALS = [\
    "Nature", "Science", "Cell", \
    "Advanced Materials", "Matter", "Nano Letters", "ACS Nano",\
    "Advanced Functional Materials", "Bioactive Materials", "Biomaterials",\
    "Small", "Chemical Engineering Journal", "Acta Biomaterialia",\
    "Advanced Healthcare Materials", "Theranostics", "Bone Research",\
    "Journal of Controlled Release", "Carbon",\
    "Carbohydrate Polymers", "Journal of Materials Chemistry A",\
    "Journal of Materials Chemistry B", "ACS Applied Materials & Interfaces",\
    "Journal of Molecular Cell Biology", "Nanoscale", \
    "Composites Part B", "Advanced Science"\
]\
\
# \uc0\u20851 \u38190 \u35789 \u35780 \u20998 \u37197 \u32622 \
KEYWORDS = \{\
    "high": ["Hydrogel", "Bone tissue engineering"], # \uc0\u26680 \u24515 \u35789  +10\u20998 \
    "medium": ["Photothermal", "Magnetothermal", "Biomaterials", "Osteogenic"], # \uc0\u30456 \u20851 \u35789  +5\u20998 \
    "negative": ["Correction", "Retraction", "Reply"] # \uc0\u25490 \u38500 \u35789 \
\}\
\
def is_target_journal(journal_name):\
    # \uc0\u27169 \u31946 \u21305 \u37197 \u65292 \u19981 \u21306 \u20998 \u22823 \u23567 \u20889 \
    j_lower = journal_name.lower()\
    for t in TARGET_JOURNALS:\
        if t.lower() in j_lower:\
            return True\
    return False\
\
def calculate_score(text):\
    text = text.lower()\
    score = 0\
    # \uc0\u25490 \u38500 \u26657 \u27491 \u12289 \u25764 \u31295 \
    for nw in KEYWORDS["negative"]:\
        if nw.lower() in text: return -100\
        \
    for kw in KEYWORDS["high"]:\
        if kw.lower() in text: score += 10\
    for kw in KEYWORDS["medium"]:\
        if kw.lower() in text: score += 5\
        \
    # \uc0\u29305 \u21035 \u21152 \u20998 \u65306 \u26082 \u26377 \u27700 \u20957 \u33014 \u21448 \u26377 \u39592 \
    if "hydrogel" in text and "bone" in text:\
        score += 15\
        \
    return score}