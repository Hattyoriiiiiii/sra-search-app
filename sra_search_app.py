from flask import Flask, render_template, request, send_file
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField, SelectMultipleField, SubmitField
from wtforms.validators import DataRequired
from pysradb.search import SraSearch
import pandas as pd
import re

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

def search_sra(organism, return_max, strategies, query):
    combined_df = pd.DataFrame()

    # columns_to_keep = [
    #     'study_accession', 'experiment_accession', 'experiment_title',
    #     'sample_taxon_id', 'sample_scientific_name', 'experiment_library_strategy',
    #     'experiment_library_source', 'experiment_library_selection', 'sample_accession',
    #     'sample_alias', 'experiment_instrument_model', 'pool_member_spots',
    #     'run_1_size', 'run_1_accession', 'run_1_total_spots', 'run_1_total_bases',
    #     'experiment_alias', 'experiment_library_construction_protocol',
    #     'experiment_link_1_type', 'experiment_link_1_value_1', 'experiment_link_1_value_2',
    #     'experiment_link_1_value_3', 'experiment_platform', 'experiment_sample_descriptor_accession',
    #     'library_layout', 'pool_external_id', 'pool_member_accession', 'pool_member_bases',
    #     'pool_member_organism', 'pool_member_sample_name', 'pool_member_sample_title',
    #     'pool_member_tax_id', 'run_1_published', 'run_1_srafile_1_date', 'run_1_srafile_1_filename',
    #     'run_1_srafile_2_filename', 'run_1_srafile_2_md5', 'run_1_srafile_2_semantic_name',
    #     'run_1_srafile_3_filename', 'run_1_srafile_3_md5', 'run_1_srafile_3_semantic_name',
    #     'sample_attributes_1_tag', 'sample_attributes_1_value', 'sample_attributes_2_tag',
    #     'sample_attributes_2_value', 'sample_attributes_3_tag', 'sample_attributes_3_value',
    #     'sample_attributes_4_tag', 'sample_attributes_4_value', 'sample_attributes_5_tag',
    #     'sample_attributes_5_value', 'sample_attributes_6_tag', 'sample_attributes_6_value',
    #     'sample_attributes_7_tag', 'sample_attributes_7_value', 'sample_external_id_1',
    #     'sample_external_id_1_namespace', 'sample_link_1_type', 'sample_link_1_value_1',
    #     'sample_link_1_value_2', 'sample_link_1_value_3', 'sample_title', 'study_alias',
    #     'study_attributes_1_tag', 'study_attributes_1_value', 'study_center_name',
    #     'study_center_project_name', 'study_external_id_1', 'study_external_id_1_namespace',
    #     'study_link_1_type', 'study_link_1_value_1', 'study_link_1_value_2',
    #     'study_study_abstract', 'study_study_title', 'study_study_type_existing_study_type',
    #     'submission_accession', 'submission_alias', 'run_1_srafile_4_filename', 'run_1_srafile_4_md5'
    # ]
    columns_to_keep = [
        'study_accession', 'experiment_accession', 'sample_scientific_name',
        'experiment_library_strategy', 'experiment_library_source', 'pubmed_id',
        'run_1_srafile_1_date', 'experiment_library_construction_protocol',
        'pool_member_sample_title', 'study_study_abstract'
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

    # ChIP-SeqとRNA-Seqの両方を含むPubMed IDを特定
    chip_seq_ids = set(combined_df[combined_df['experiment_library_strategy'] == 'ChIP-Seq']['pubmed_id'])
    rna_seq_ids = set(combined_df[combined_df['experiment_library_strategy'] == 'RNA-Seq']['pubmed_id'])
    valid_ids = chip_seq_ids.intersection(rna_seq_ids)

    # 有効なPubMed IDに一致する行を保持
    filtered_df = combined_df[combined_df['pubmed_id'].isin(valid_ids)]

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
