import streamlit as st
import pandas as pd
import json
import os

st.set_page_config(page_title="BioMed Daily", page_icon="🧪", layout="wide")

st.title("🧪 每日生物医学文献推送")
st.caption("Focus: Hydrogel | Bone Tissue | Biomaterials (JMCB to CNS)")

# 读取数据
if os.path.exists("data.json"):
    with open("data.json", "r", encoding='utf-8') as f:
        data = json.load(f)
    df = pd.DataFrame(data)
else:
    df = pd.DataFrame()

if not df.empty:
    # 侧边栏
    min_score = st.sidebar.slider("最低评分筛选", 0, 50, 10)
    search = st.sidebar.text_input("关键词搜索")
    
    # 筛选
    df = df[df['score'] >= min_score]
    if search:
        df = df[df['title'].str.contains(search, case=False) | df['abstract'].str.contains(search, case=False)]
    
    st.sidebar.markdown(f"共显示 {len(df)} 篇文献")

    # 展示
    for i, row in df.iterrows():
        with st.container():
            st.markdown(f"### [{row['title']}]({row['url']})")
            st.markdown(f"**{row['journal']}** | {row['date']} | Score: `{row['score']}`")
            
            try:
                ai_data = json.loads(row['analysis'])
            except:
                ai_data = {}
            
            col1, col2 = st.columns(2)
            with col1:
                st.success(f"**中文摘要**: {ai_data.get('chinese_summary', '无')}")
            with col2:
                st.info(f"**创新点**: {ai_data.get('innovation', '无')}")
                st.warning(f"**不足**: {ai_data.get('flaw', '无')}")
            
            with st.expander("查看原文摘要"):
                st.write(row['abstract'])
            st.divider()
else:
    st.info("暂无数据，请等待 GitHub Action 首次运行。")