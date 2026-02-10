import streamlit as st
import pandas as pd
import json
import os

# 页面配置
st.set_page_config(page_title="BioMed Daily", page_icon="🧬", layout="wide")

st.title("🧬 每日生物医学文献推送")
st.caption("Focus: Hydrogel | Bone Tissue | Biomaterials (JMCB to CNS)")

# --- 核心修改：读取 JSON 文件 ---
def load_data():
    if os.path.exists("data.json"):
        try:
            with open("data.json", "r", encoding='utf-8') as f:
                data = json.load(f)
            # 转换为 DataFrame
            return pd.DataFrame(data)
        except Exception as e:
            st.error(f"读取数据出错: {e}")
            return pd.DataFrame()
    else:
        return pd.DataFrame()

df = load_data()

# --- 展示逻辑 ---
if not df.empty:
    # 确保 score 列是数字类型，防止排序出错
    if 'score' in df.columns:
        df['score'] = pd.to_numeric(df['score'], errors='coerce').fillna(0)
    
    # 侧边栏筛选
    with st.sidebar:
        st.header("🔍 筛选")
        min_score = st.slider("最低评分", 0, 50, 10)
        search_text = st.text_input("搜索标题/摘要")
        st.markdown(f"**共找到 {len(df)} 篇文献**")

    # 应用筛选
    filtered_df = df[df['score'] >= min_score]
    
    if search_text:
        filtered_df = filtered_df[
            filtered_df['title'].str.contains(search_text, case=False) | 
            filtered_df['abstract'].str.contains(search_text, case=False)
        ]

    # 按日期和分数降序排列
    filtered_df = filtered_df.sort_values(by=['date', 'score'], ascending=[False, False])

    # 循环展示每一篇
    for i, row in filtered_df.iterrows():
        with st.container():
            # 标题链接
            st.markdown(f"### [{row['title']}]({row['url']})")
            
            # 元数据行
            st.markdown(
                f"**📅 {row['date']}** | 📖 *{row['journal']}* | 🔥 Score: `{row['score']}`"
            )

            # 解析 AI 分析结果
            try:
                # 如果 analysis 已经是字典就直接用，如果是字符串就转一下
                if isinstance(row.get('analysis'), dict):
                    ai_data = row['analysis']
                else:
                    ai_data = json.loads(row.get('analysis', '{}'))
            except:
                ai_data = {}

            # 两列布局展示 AI 结果
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.success(f"**📝 中文摘要**: {ai_data.get('chinese_summary', '暂无')}")
            
            with col2:
                st.info(f"**💡 创新点**: {ai_data.get('innovation', '暂无')}")
                st.warning(f"**⚠️ 不足/局限**: {ai_data.get('flaw', '暂无')}")
            
            # 原文摘要折叠
            with st.expander("查看英文原摘要"):
                st.write(row.get('abstract', 'No abstract available.'))
            
            st.markdown("---")

else:
    st.info("📭 暂无数据。请检查 GitHub Actions 是否成功运行并生成了 data.json 文件。")
