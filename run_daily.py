{\rtf1\ansi\ansicpg936\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 import os\
import json\
import datetime\
import requests\
from Bio import Entrez\
from openai import OpenAI\
from config import KEYWORDS, is_target_journal, calculate_score\
\
# \uc0\u37197 \u32622 \
DATA_FILE = "data.json"\
DEEPSEEK_KEY = os.environ.get("DEEPSEEK_API_KEY") # \uc0\u20174 \u29615 \u22659 \u21464 \u37327 \u33719 \u21462 Key\
Entrez.email = os.environ.get("USER_EMAIL") # \uc0\u20174 \u29615 \u22659 \u21464 \u37327 \u33719 \u21462 \u37038 \u31665 \
\
def fetch_pubmed_today():\
    print("\uc0\u27491 \u22312 \u25628 \u32034  PubMed...")\
    # \uc0\u33719 \u21462 \u26368 \u36817 2\u22825 \u30340 \u25991 \u29486 \
    date_query = '("2024/01/01"[Date - Publication] : "3000"[Date - Publication])'\
    # \uc0\u20851 \u38190 \u35789 \u26500 \u24314 \
    term_query = '("Hydrogels"[Mesh] OR "Bone Tissue Engineering"[Mesh] OR "Biomaterials"[Mesh])'\
    query = f"\{term_query\} AND \{date_query\}"\
    \
    try:\
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="date")\
        record = Entrez.read(handle)\
        id_list = record["IdList"]\
        if not id_list: return []\
        \
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")\
        papers = Entrez.read(handle)\
        \
        candidates = []\
        for paper in papers['PubmedArticle']:\
            art = paper['MedlineCitation']['Article']\
            journal = art['Journal']['Title']\
            title = art['ArticleTitle']\
            \
            # 1. \uc0\u26399 \u21002 \u21021 \u31579 \
            if not is_target_journal(journal):\
                continue\
                \
            abstract = art['Abstract']['AbstractText'][0] if 'Abstract' in art else ""\
            \
            # 2. \uc0\u20851 \u38190 \u35789 \u35780 \u20998 \
            score = calculate_score(title + " " + abstract)\
            if score < 10: # \uc0\u20998 \u25968 \u22826 \u20302 \u19981 \u35201 \
                continue\
                \
            doi = ""\
            for x in paper['PubmedData']['ArticleIdList']:\
                if x.attributes.get('IdType') == 'doi':\
                    doi = str(x)\
            \
            candidates.append(\{\
                "id": doi if doi else title, # \uc0\u29992 DOI\u20570 \u21807 \u19968 ID\
                "title": title,\
                "journal": journal,\
                "authors": ", ".join([a.get('LastName','') for a in art.get('AuthorList',[])][:3]) + " et al.",\
                "abstract": abstract,\
                "score": score,\
                "url": f"https://pubmed.ncbi.nlm.nih.gov/\{paper['MedlineCitation']['PMID']\}/",\
                "date": datetime.datetime.now().strftime("%Y-%m-%d"),\
                "is_new": True\
            \})\
        return candidates\
    except Exception as e:\
        print(f"Error fetching: \{e\}")\
        return []\
\
def analyze_article(text):\
    if not DEEPSEEK_KEY: return "Error: No API Key"\
    \
    client = OpenAI(api_key=DEEPSEEK_KEY, base_url="https://api.deepseek.com")\
    prompt = f"""\
    \uc0\u35831 \u20998 \u26512 \u36825 \u31687 \u29983 \u29289 \u21307 \u23398 \u25991 \u29486 \u25688 \u35201 \u12290 \
    \uc0\u25688 \u35201 : \{text\}\
    \
    \uc0\u35831 \u36820 \u22238 \u19968 \u20010 JSON\u26684 \u24335 \u23383 \u31526 \u20018 (\u19981 \u35201 Markdown)\u65292 \u21253 \u21547 \u20197 \u19979 \u23383 \u27573 \u65306 \
    chinese_summary: \uc0\u20013 \u25991 \u25688 \u35201 (\u31616 \u32451 )\
    innovation: \uc0\u21019 \u26032 \u28857 (1\u21477 \u35805 )\
    flaw: \uc0\u28508 \u22312 \u19981 \u36275 \u25110 \u23616 \u38480 (1\u21477 \u35805 )\
    """\
    try:\
        response = client.chat.completions.create(\
            model="deepseek-chat",\
            messages=[\{"role": "user", "content": prompt\}],\
            response_format=\{ "type": "json_object" \}\
        )\
        return response.choices[0].message.content\
    except:\
        return "\{\}"\
\
def main():\
    # 1. \uc0\u35835 \u21462 \u26087 \u25968 \u25454 \
    if os.path.exists(DATA_FILE):\
        with open(DATA_FILE, 'r', encoding='utf-8') as f:\
            old_data = json.load(f)\
    else:\
        old_data = []\
    \
    old_ids = [d['id'] for d in old_data]\
    \
    # 2. \uc0\u33719 \u21462 \u26032 \u25968 \u25454 \
    candidates = fetch_pubmed_today()\
    new_data = []\
    \
    for item in candidates:\
        if item['id'] not in old_ids:\
            print(f"DeepSeek \uc0\u20998 \u26512 \u20013 : \{item['title'][:20]\}...")\
            analysis = analyze_article(item['abstract'])\
            item['analysis'] = analysis\
            new_data.append(item)\
    \
    # 3. \uc0\u20445 \u23384 \
    if new_data:\
        print(f"\uc0\u26032 \u22686  \{len(new_data)\} \u31687 \u25991 \u29486 ")\
        final_data = new_data + old_data\
        with open(DATA_FILE, 'w', encoding='utf-8') as f:\
            json.dump(final_data, f, ensure_ascii=False, indent=2)\
    else:\
        print("\uc0\u26080 \u26032 \u25991 \u29486 \u26356 \u26032 ")\
\
if __name__ == "__main__":\
    main()}