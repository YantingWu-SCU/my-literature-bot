# 核心期刊列表 (JMCB - CNS 级别)
# 只要期刊名包含以下字符串，就会被保留
TARGET_JOURNALS = [
    "Nature", "Science", "Cell", 
    "Advanced Materials", "Matter", "Nano Letters", "ACS Nano",
    "Advanced Functional Materials", "Bioactive Materials", "Biomaterials",
    "Small", "Chemical Engineering Journal", "Acta Biomaterialia",
    "Advanced Healthcare Materials", "Theranostics", "Bone Research",
    "Journal of Controlled Release", "Carbon",
    "Carbohydrate Polymers", "Journal of Materials Chemistry A",
    "Journal of Materials Chemistry B", "ACS Applied Materials & Interfaces",
    "Journal of Molecular Cell Biology", "Nanoscale", 
    "Composites Part B", "Advanced Science"
]

# 关键词评分配置
KEYWORDS = {
    "high": ["Hydrogel", "Bone tissue engineering"], # 核心词 +10分
    "medium": ["Photothermal", "Magnetothermal", "Biomaterials", "Osteogenic"], # 相关词 +5分
    "negative": ["Correction", "Retraction", "Reply"] # 排除词
}

def is_target_journal(journal_name):
    # 模糊匹配，不区分大小写
    j_lower = journal_name.lower()
    for t in TARGET_JOURNALS:
        if t.lower() in j_lower:
            return True
    return False

def calculate_score(text):
    text = text.lower()
    score = 0
    # 排除校正、撤稿
    for nw in KEYWORDS["negative"]:
        if nw.lower() in text: return -100
        
    for kw in KEYWORDS["high"]:
        if kw.lower() in text: score += 10
    for kw in KEYWORDS["medium"]:
        if kw.lower() in text: score += 5
        
    # 特别加分：既有水凝胶又有骨
    if "hydrogel" in text and "bone" in text:
        score += 15
        
    return score
