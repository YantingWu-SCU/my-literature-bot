import os
import json
import datetime
import requests
from Bio import Entrez
from openai import OpenAI
from config import KEYWORDS, is_target_journal, calculate_score

# 配置
DATA_FILE = "data.json"
DEEPSEEK_KEY = os.environ.get("DEEPSEEK_API_KEY") # 从环境变量获取Key
Entrez.email = os.environ.get("USER_EMAIL") # 从环境变量获取邮箱

def fetch_pubmed_today():
    print("正在搜索 PubMed...")
    # 获取最近2天的文献
    date_query = '("2024/01/01"[Date - Publication] : "3000"[Date - Publication])'
    # 关键词构建
    term_query = '("Hydrogels"[Mesh] OR "Bone Tissue Engineering"[Mesh] OR "Biomaterials"[Mesh])'
    query = f"{term_query} AND {date_query}"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="date")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        papers = Entrez.read(handle)
        
        candidates = []
        for paper in papers['PubmedArticle']:
            art = paper['MedlineCitation']['Article']
            journal = art['Journal']['Title']
            title = art['ArticleTitle']
            
            # 1. 期刊初筛
            if not is_target_journal(journal):
                continue
                
            abstract = art['Abstract']['AbstractText'][0] if 'Abstract' in art else ""
            
            # 2. 关键词评分
            score = calculate_score(title + " " + abstract)
            if score < 10: # 分数太低不要
                continue
                
            doi = ""
            for x in paper['PubmedData']['ArticleIdList']:
                if x.attributes.get('IdType') == 'doi':
                    doi = str(x)
            
            candidates.append({
                "id": doi if doi else title, # 用DOI做唯一ID
                "title": title,
                "journal": journal,
                "authors": ", ".join([a.get('LastName','') for a in art.get('AuthorList',[])][:3]) + " et al.",
                "abstract": abstract,
                "score": score,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{paper['MedlineCitation']['PMID']}/",
                "date": datetime.datetime.now().strftime("%Y-%m-%d"),
                "is_new": True
            })
        return candidates
    except Exception as e:
        print(f"Error fetching: {e}")
        return []

def analyze_article(text):
    if not DEEPSEEK_KEY: return "Error: No API Key"
    
    client = OpenAI(api_key=DEEPSEEK_KEY, base_url="https://api.deepseek.com")
    prompt = f"""
    请分析这篇生物医学文献摘要。
    摘要: {text}
    
    请返回一个JSON格式字符串(不要Markdown)，包含以下字段：
    chinese_summary: 中文摘要(简练)
    innovation: 创新点(1句话)
    flaw: 潜在不足或局限(1句话)
    """
    try:
        response = client.chat.completions.create(
            model="deepseek-chat",
            messages=[{"role": "user", "content": prompt}],
            response_format={ "type": "json_object" }
        )
        return response.choices[0].message.content
    except:
        return "{}"

def main():
    # 1. 读取旧数据
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, 'r', encoding='utf-8') as f:
            old_data = json.load(f)
    else:
        old_data = []
    
    old_ids = [d['id'] for d in old_data]
    
    # 2. 获取新数据
    candidates = fetch_pubmed_today()
    new_data = []
    
    for item in candidates:
        if item['id'] not in old_ids:
            print(f"DeepSeek 分析中: {item['title'][:20]}...")
            analysis = analyze_article(item['abstract'])
            item['analysis'] = analysis
            new_data.append(item)
    
    # 3. 保存
    if new_data:
        print(f"新增 {len(new_data)} 篇文献")
        final_data = new_data + old_data
        with open(DATA_FILE, 'w', encoding='utf-8') as f:
            json.dump(final_data, f, ensure_ascii=False, indent=2)
    else:
        print("无新文献更新")

if __name__ == "__main__":
    main()
