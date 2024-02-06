from flask import Flask, render_template, request, send_file
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField, SelectMultipleField, SubmitField
from wtforms.validators import DataRequired
from pysradb.search import SraSearch
import pandas as pd
import re
from Bio import Entrez

Entrez.email = "www.tatsuya92@gmail.com"

app = Flask(__name__, template_folder='templates')
app.config['SECRET_KEY'] = 'your_secret_key'
download_file_name = ""

class SearchForm(FlaskForm):
    organism = StringField('Organism', validators=[DataRequired()], default='human')
    return_max = IntegerField('Return Max', validators=[DataRequired()], default=10000)
    strategy = SelectMultipleField('Strategy', 
        choices=[(strategy, strategy) for strategy in ['AMPLICON', 'ATAC-seq', 'Bisulfite-Seq', 
                 'CLONE', 'CLONEEND', 'CTS', 'ChIA-PET', 'ChIP-Seq', 'DNase-Hypersensitivity', 
                 'EST', 'FAIRE-seq', 'FINISHING', 'FL-cDNA', 'Hi-C', 'MBD-Seq', 
                 'MNase-Seq', 'MRE-Seq', 'MeDIP-Seq', 'OTHER', 'POOLCLONE', 'RAD-Seq', 
                 'RIP-Seq', 'RNA-Seq', 'SELEX', 'Synthetic-Long-Read', 'Targeted-Capture', 
                 'Tethered Chromatin Conformation Capture', 'Tn-Seq', 'VALIDATION', 'WCS', 
                 'WGA', 'WGS', 'WXS', 'miRNA-Seq', 'ncRNA-Seq', 'ssRNA-seq', 'GBS']],
        validators=[DataRequired()], default=['RNA-Seq', 'ChIP-Seq'])
    query = StringField('Query')
    submit = SubmitField('Search')

def format_link(id):
    return f'<a href="https://pubmed.ncbi.nlm.nih.gov/{id}">{id}</a>'

def fetch_pubmed_details(pubmed_id):
    handle = Entrez.efetch(db="pubmed", id=str(pubmed_id), retmode="xml")
    records = Entrez.read(handle)
    details = {
        'title': None,
        'journal': None,
        'publication_date': None
    }
    try:
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        details['title'] = article.get('ArticleTitle', None)
        if 'Journal' in article:
            details['journal'] = article['Journal'].get('Title', None)
        if 'ArticleDate' in article and article['ArticleDate']:
            details['publication_date'] = f"{article['ArticleDate'][0]['Year']}-{article['ArticleDate'][0]['Month']}-{article['ArticleDate'][0].get('Day', '01')}"
        elif 'PubDate' in article and article['PubDate']:
            pub_date = article['PubDate']
            year = pub_date.get('Year', 'Unknown')
            month = pub_date.get('Month', '01')
            day = pub_date.get('Day', '01')
            details['publication_date'] = f"{year}-{month}-{day}"
    except IndexError as e:
        print(f"Error fetching details for PubMed ID {pubmed_id}: {e}")
    return details


def update_dataframe_with_pubmed_info(df):
    titles, journals, publication_dates = [], [], []
    for pubmed_id in df['pubmed_id']:
        if pd.notnull(pubmed_id):
            details = fetch_pubmed_details(pubmed_id)
            titles.append(details['title'])
            journals.append(details['journal'])
            publication_dates.append(details['publication_date'])
        else:
            titles.append(None)
            journals.append(None)
            publication_dates.append(None)
    
    df['title'] = titles
    df['journal'] = journals
    df['publication_date'] = publication_dates
    return df

def add_h3k27ac_flag(df):
    # H3K27acをターゲットにしているかどうかを示すフラグを追加
    df['h3k27ac_flag'] = df.apply(lambda row: 1 if 'ChIP-Seq' in row['experiment_library_strategy'] and 'H3K27ac' in row.get('pool_member_sample_title', '') else 0, axis=1)
    return df


def search_sra(organism, return_max, strategies, query):
    combined_df = pd.DataFrame()

    # columns_to_keep = [
    #     'study_accession', 'experiment_accession', 'sample_scientific_name',
    #     'experiment_library_strategy', 'experiment_library_source', 'pubmed_id',
    #     'run_1_srafile_1_date', 'experiment_library_construction_protocol',
    #     'pool_member_sample_title', 'study_study_abstract'
    # ]
    columns_to_keep = [
        'study_external_id_1', 'run_1_accession', 
        'pool_member_sample_title', 'study_study_abstract', 'study_study_title',
        'pool_member_organism', 'experiment_instrument_model', 'experiment_library_strategy',
        'experiment_library_source', 'experiment_library_selection', 'library_layout',
        'experiment_library_construction_protocol', 'study_alias', 'run_1_published', 'pubmed_id'
    ]

    for strategy in strategies:
        instance = SraSearch(
            verbosity=3, 
            return_max=return_max, 
            organism=organism, 
            strategy=strategy, 
            query=query
        )
        instance.search()
        df = instance.get_df()
        combined_df = pd.concat([combined_df, df])

    # 各行に対してPubMed IDを抽出
    def extract_pubmed_ids(row):
        for col in combined_df.columns:
            if re.match(r'study_link_\d+_value_1', col):
                link_number = col.split('_')[2]
                link_value_2_col = f'study_link_{link_number}_value_2'
                if pd.notna(row[col]) and row[col] == 'DB: pubmed' and link_value_2_col in combined_df.columns:
                    pubmed_id = row[link_value_2_col].replace('ID: ', '')
                    return int(pubmed_id) if pubmed_id.isdigit() else None
        return None

    combined_df['pubmed_id'] = combined_df.apply(extract_pubmed_ids, axis=1)
    combined_df = combined_df[columns_to_keep]
    combined_df = add_h3k27ac_flag(combined_df)

    # ChIP-SeqとRNA-Seqの両方を含むPubMed IDを特定
    chip_seq_ids = set(combined_df[combined_df['experiment_library_strategy'] == 'ChIP-Seq']['pubmed_id'])
    rna_seq_ids = set(combined_df[combined_df['experiment_library_strategy'] == 'RNA-Seq']['pubmed_id'])
    valid_ids = chip_seq_ids.intersection(rna_seq_ids)

    # 有効なPubMed IDに一致する行を保持
    filtered_df = combined_df[combined_df['pubmed_id'].isin(valid_ids)]

    # PubMedからのデータを取得
    filtered_df = update_dataframe_with_pubmed_info(filtered_df)

    # 各PubMed IDについて選択されたstrategyの数を計算
    strategy_count = filtered_df.groupby('pubmed_id')['experiment_library_strategy'].value_counts().unstack().fillna(0)
    strategy_count = strategy_count.astype(int)

    # リンクを追加
    strategy_count.index = strategy_count.index.astype(int).map(format_link)

    return filtered_df, strategy_count

@app.route('/', methods=['GET', 'POST'])
def index():
    global download_file_name  # グローバル変数を使用

    form = SearchForm()
    if form.validate_on_submit():
        organism = form.organism.data
        return_max = form.return_max.data
        strategies = form.strategy.data
        query = form.query.data

        filtered_df, strategy_count = search_sra(organism, return_max, strategies, query)

        # クエリをファイル名に変換して保存
        download_file_name = query.replace(" ", "-").lower() + ".tsv"
        filtered_df.to_csv(download_file_name, sep="\t", index=False)

        return render_template('results.html', 
                           results_table=filtered_df.to_html(classes='data'), 
                           strategy_counts=strategy_count.to_html(classes='data', escape=False),
                           titles=filtered_df.columns.values)
    return render_template('index.html', form=form)

@app.route('/download')
def download():
    global download_file_name  # グローバル変数を使用
    return send_file(download_file_name, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
