import os
import json
import datetime
import time
from Bio import Entrez
from openai import OpenAI
from config import KEYWORDS, is_target_journal, calculate_score

DATA_FILE = "data.json"
# 获取环境变量
DEEPSEEK_KEY = os.environ.get("DEEPSEEK_API_KEY")
Entrez.email = os.environ.get("USER_EMAIL")

def analyze_article(text):
    # 1. 检查 Key 是否存在
    if not DEEPSEEK_KEY:
        print("Error: 找不到 DEEPSEEK_API_KEY，请检查 GitHub Secrets！")
        return json.dumps({"chinese_summary": "API Key 未配置", "innovation": "无法分析", "flaw": "请在 GitHub Settings 添加 Key"})

    client = OpenAI(api_key=DEEPSEEK_KEY, base_url="https://api.deepseek.com")
    
    prompt = f"""
    你是一个生物医学专家。请分析以下文献摘要：
    {text[:2000]} 
    
    请务必返回合法的 JSON 格式，包含以下三个字段：
    chinese_summary: 中文摘要（简练）
    innovation: 创新点（1句话）
    flaw: 不足或局限（1句话）
    """
    
    try:
        response = client.chat.completions.create(
            model="deepseek-chat",
            messages=[
                {"role": "system", "content": "You are a helpful assistant. You must response in JSON format."},
                {"role": "user", "content": prompt}
            ],
            response_format={ "type": "json_object" },
            temperature=0.1
        )
        return response.choices[0].message.content
    except Exception as e:
        # 重点：打印具体报错信息到日志
        print(f"DeepSeek API Error: {e}")
        # 返回一个包含错误信息的 JSON，方便在网页显示
        return json.dumps({
            "chinese_summary": f"分析出错: {str(e)[:20]}...", 
            "innovation": "API 调用失败", 
            "flaw": "请检查 GitHub Logs"
        })

def fetch_pubmed_today():
    print("正在搜索 PubMed...")
    # 稍微放宽一点日期，确保能抓到数据测试
    date_query = '("2024/01/01"[Date - Publication] : "3000"[Date - Publication])'
    term_query = '("Hydrogels"[Mesh] OR "Bone Tissue Engineering"[Mesh] OR "Biomaterials"[Mesh])'
    query = f"{term_query} AND {date_query}"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10, sort="date")
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
            
            # 期刊筛选
            if not is_target_journal(journal):
                continue
            
            abstract = art['Abstract']['AbstractText'][0] if 'Abstract' in art else "No Abstract"
            
            # 评分筛选
            score = calculate_score(title + " " + abstract)
            if score < 5: continue 
                
            doi = ""
            for x in paper['PubmedData']['ArticleIdList']:
                if x.attributes.get('IdType') == 'doi':
                    doi = str(x)
            
            candidates.append({
                "id": doi if doi else title,
                "title": title,
                "journal": journal,
                "authors": ", ".join([a.get('LastName','') for a in art.get('AuthorList',[])][:3]),
                "abstract": abstract,
                "score": score,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{paper['MedlineCitation']['PMID']}/",
                "date": datetime.datetime.now().strftime("%Y-%m-%d"),
            })
        return candidates
    except Exception as e:
        print(f"Fetch Error: {e}")
        return []

def main():
    if os.path.exists(DATA_FILE):
        try:
            with open(DATA_FILE, 'r', encoding='utf-8') as f:
                old_data = json.load(f)
        except:
            old_data = []
    else:
        old_data = []
    
    old_ids = [d['id'] for d in old_data]
    
    # 获取新数据
    candidates = fetch_pubmed_today()
    new_data = []
    
    for item in candidates:
        if item['id'] not in old_ids:
            print(f"正在分析 AI: {item['title'][:30]}...")
            analysis = analyze_article(item['abstract'])
            item['analysis'] = analysis
            new_data.append(item)
            # 避免请求过快
            time.sleep(2)
    
    if new_data:
        print(f"存入 {len(new_data)} 篇新文献")
        # 新数据放在最前面
        final_data = new_data + old_data
        with open(DATA_FILE, 'w', encoding='utf-8') as f:
            json.dump(final_data, f, ensure_ascii=False, indent=2)
    else:
        print("没有符合条件的新文献")

if __name__ == "__main__":
    main()
