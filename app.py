# app.py
import streamlit as st
import sqlite3
import pandas as pd
import json

# 设置页面风格
st.set_page_config(page_title="BioMed Research Feed", layout="wide", page_icon="🧬")

# 自定义 CSS 以接近 ObservableHQ 风格
st.markdown("""
<style>
    .reportview-container { background: #ffffff; }
    .main-card { 
        padding: 20px; 
        border-radius: 10px; 
        border: 1px solid #e0e0e0; 
        margin-bottom: 20px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .tag {
        display: inline-block;
        padding: 2px 8px;
        border-radius: 4px;
        font-size: 0.8em;
        margin-right: 5px;
    }
    .tag-review { background-color: #e3f2fd; color: #1565c0; }
    .tag-research { background-color: #e8f5e9; color: #2e7d32; }
    .tag-score { background-color: #fff3e0; color: #ef6c00; font-weight: bold; }
    h3 { color: #2c3e50; font-family: 'Helvetica Neue', sans-serif; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 生物医学文献智能推送系统")
st.caption("Focus: Hydrogels, Bone Tissue Engineering, Biomaterials | Powered by DeepSeek")

# 读取数据
@st.cache_data(ttl=3600) # 缓存数据减少IO
def load_data():
    conn = sqlite3.connect("literature.db")
    df = pd.read_sql_query("SELECT * FROM articles ORDER BY fetch_date DESC, keyword_score DESC", conn)
    conn.close()
    return df

df = load_data()

# --- 侧边栏筛选 ---
st.sidebar.header("🔍 筛选条件")

# 1. 文章类型筛选
types = st.sidebar.multiselect("文章类型", options=df['article_type'].unique(), default=df['article_type'].unique())

# 2. 评分筛选
min_score = st.sidebar.slider("最低关键词评分", int(df['keyword_score'].min()), int(df['keyword_score'].max()), 5)

# 3. 来源筛选
sources = st.sidebar.multiselect("来源", options=df['source'].unique(), default=df['source'].unique())

# 4. 搜索框
search_text = st.sidebar.text_input("搜索标题/作者")

# 应用筛选
filtered_df = df[
    (df['article_type'].isin(types)) & 
    (df['keyword_score'] >= min_score) &
    (df['source'].isin(sources))
]

if search_text:
    filtered_df = filtered_df[filtered_df['title'].str.contains(search_text, case=False)]

st.sidebar.markdown(f"**展示 {len(filtered_df)} / {len(df)} 篇文献**")

# --- 主展示区 ---

for index, row in filtered_df.iterrows():
    # 解析 DeepSeek 返回的 JSON
    try:
        analysis = json.loads(row['deepseek_analysis'])
    except:
        analysis = {"chinese_translation": "解析失败", "main_findings": "无", "innovation": "无", "flaws": "无", "future_directions": "无"}

    # 构建卡片 UI
    with st.container():
        st.markdown(f"""
        <div class="main-card">
            <h3><a href="{row['url']}" target="_blank" style="text-decoration:none; color:#2c3e50;">{row['title']}</a></h3>
            <p style="color: #666; font-size: 0.9em;">
                <b>{row['journal']}</b> (IF: {row['if_score']}) | {row['year']} | {row['authors']}
            </p>
            <div>
                <span class="tag tag-{'review' if row['article_type']=='Review' else 'research'}">{row['article_type']}</span>
                <span class="tag tag-score">🔥 Score: {row['keyword_score']}</span>
                <span class="tag">{row['source']}</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        with st.expander("🤖 DeepSeek 深度分析 & 中文摘要"):
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown("#### 📝 中文摘要")
                st.write(analysis.get("chinese_translation", "暂无翻译"))
                
                st.markdown("#### 💡 创新点")
                st.info(analysis.get("innovation", "暂无"))
                
            with col2:
                st.markdown("#### 🔍 主要发现")
                st.success(analysis.get("main_findings", "暂无"))
                
                st.markdown("#### ⚠️ 缺陷/局限性")
                st.warning(analysis.get("flaws", "暂无"))
                
                st.markdown("#### 🚀 未来方向")
                st.write(analysis.get("future_directions", "暂无"))
            
            st.markdown("---")
            st.markdown("**原始摘要:**")
            st.caption(row['abstract'])

    st.write("") # Spacer
