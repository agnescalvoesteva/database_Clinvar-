# The script of python that you have to run to create the website with the database
from flask import Flask, render_template, request, redirect, url_for, jsonify
from markupsafe import Markup
import sqlite3
import pandas as pd
import numpy as np
import requests
import time
import logging

# Monkey patch para flask-paginate
import flask_paginate
flask_paginate.Markup = Markup

from flask_paginate import Pagination, get_page_parameter

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

app = Flask(__name__)
DB_FILE = "clinvar_db.sqlite"

# -------------------- CLINVAR SEARCH FUNCTIONS --------------------
def search_variants_by_condition(condition, max_results=50):
    """Search ClinVar for variants by condition"""
    print(f"Looking up variants for: {condition}")
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": condition,
        "retmax": max_results,
        "retmode": "json"
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        data = response.json()
        variant_ids = data.get("esearchresult", {}).get("idlist", [])
        
        if not variant_ids:
            print(f"No variants found for: {condition}")
            return pd.DataFrame()
        
        print(f"Found {len(variant_ids)} variants")
        return get_variant_details(variant_ids)
        
    except Exception as e:
        print(f"Search failed: {e}")
        return pd.DataFrame()

# Función corregida para consultar ClinVar
def query_clinvar(condition, filters=None):
    """Query ClinVar API with optional filters"""
    print(f"Querying ClinVar for: {condition} with filters: {filters}")
    
    # Construir término de búsqueda
    search_term = condition
    
    if filters and filters.get('variation_type'):
        var_type = filters['variation_type']
        # Mapear a términos de búsqueda de ClinVar
        type_map = {
            'Deletion': 'deletion',
            'Duplication': 'duplication',
            'Indel': 'indel',
            'Insertion': 'insertion',
            'Single nucleotide': 'single nucleotide variant',
            'Inversion': 'inversion',
            'Translocation': 'translocation',
            'Tandem duplication': 'tandem duplication',
            'Copy number loss': 'copy number loss',
            'Complex': 'complex'
        }
        if var_type in type_map:
            search_term = f"{condition} AND {type_map[var_type]}[variant type]"
    
    # Llamar a la función existente search_variants_by_condition
    max_results = 50  # Puedes ajustar esto según sea necesario
    variants_df = search_variants_by_condition(search_term, max_results)
    
    # Procesar datos para la respuesta
    if variants_df.empty:
        return {
            'associated_genes': [],
            'total_variants': 0,
            'unique_genes': 0
        }
    
    # Contar genes únicos
    unique_genes = variants_df['gene_symbol'].nunique()
    
    # Preparar lista de genes con conteos
    gene_counts = variants_df['gene_symbol'].value_counts().reset_index()
    gene_counts.columns = ['gene', 'variant_count']
    
    associated_genes = []
    for _, row in gene_counts.iterrows():
        if row['gene'] and str(row['gene']).strip():
            gene_info = {
                'gene': row['gene'],
                'variant_count': int(row['variant_count'])
            }
            associated_genes.append(gene_info)
    
    return {
        'associated_genes': associated_genes,
        'total_variants': len(variants_df),
        'unique_genes': unique_genes
    }

def get_variant_details(variant_ids):
    """Get detailed information for a list of variant IDs"""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    all_variants = []
    
    for i in range(0, len(variant_ids), 50):
        batch = variant_ids[i:i+50]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json"
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            data = response.json().get("result", {})
            
            for variant_id in batch:
                if variant_id in data:
                    variant = extract_variant_info(data[variant_id])
                    all_variants.append(variant)
            
        except Exception as e:
            print(f"Batch error: {e}")
        
        time.sleep(0.5)
    
    return pd.DataFrame(all_variants)

def extract_variant_info(variant_data):
    """Extract relevant information from variant data"""
    try:
        variant_id = variant_data.get('uid', '')
        title = variant_data.get('title', '')
        
        gene = ""
        if 'genes' in variant_data:
            genes = variant_data.get('genes', [])
            if genes and isinstance(genes, list):
                gene = genes[0].get('symbol', '')
        
        significance = ""
        if 'clinical_significance' in variant_data:
            clin_sig = variant_data['clinical_significance']
            if 'description' in clin_sig:
                significance = clin_sig['description']
                
        diseases = ""
        if 'trait_set' in variant_data:
            trait_set = variant_data['trait_set']
            if 'trait' in trait_set:
                traits = trait_set['trait']
                if isinstance(traits, list):
                    names = [trait.get('trait_name', '') for trait in traits if trait.get('trait_name')]
                    diseases = ", ".join(names)
        return {
            "variant_id": variant_id,
            "title": title,
            "gene_symbol": gene,
            "clinical_significance": significance,
            "conditions": diseases,
            "last_updated": variant_data.get('last_updated', '')
        }
        
    except Exception as e:
        print(f"Error parsing variant: {e}")
        return {
            "variant_id": variant_data.get('uid', ''),
            "title": variant_data.get('title', ''),
            "gene_symbol": "",
            "clinical_significance": "",
            "conditions": "",
            "last_updated": ""
        }

def get_gene_info(gene_symbol):
    """Get gene information from Ensembl"""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            ensembl_id = data.get('id')
            description = data.get('description', '')
            
            uniprot_id = None
            xrefs = data.get('xrefs', [])
            for xref in xrefs:
                if xref.get('dbname') in ['UniProtKB/Swiss-Prot', 'UniProtKB']:
                    uniprot_id = xref.get('primary_id')
                    break
            
            return ensembl_id, uniprot_id, description
        else:
            return None, None, f"Not found (HTTP {response.status_code})"
            
    except Exception as e:
        print(f"Error getting gene info for {gene_symbol}: {e}")
        return None, None, f"Error: {str(e)}"

def create_gene_data(variants_df):
    """Create gene data from variants dataframe"""
    if variants_df.empty:
        return pd.DataFrame()
    
    unique_genes = variants_df['gene_symbol'].unique()
    print(f"Processing {len(unique_genes)} genes")
    
    gene_list = []
    for gene in unique_genes:
        if gene and gene.strip():
            ensembl_id, uniprot_id, description = get_gene_info(gene)
            variant_count = len(variants_df[variants_df['gene_symbol'] == gene])
            
            gene_list.append({
                'Gene_Symbol': gene,
                'Ensembl_Gene_ID': ensembl_id,
                'UniProtKB': uniprot_id,
                'Description': description,
                'Variant_Count': variant_count
            })
            
            time.sleep(0.3)
    
    return pd.DataFrame(gene_list)

def save_variants(variants_df, condition):
    """Save variants to database"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    for _, row in variants_df.iterrows():
        cursor.execute("""
            INSERT INTO variants (variant_id, gene, clinical_significance, diseases, last_updated, searched_condition)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            str(row['variant_id']),
            row['gene_symbol'],
            row['clinical_significance'],
            row['conditions'],
            row['last_updated'],
            condition
        ))
    
    conn.commit()
    conn.close()
    print(f"Saved {len(variants_df)} variants for '{condition}'")

def save_genes(genes_df):
    """Save genes to database"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    for _, row in genes_df.iterrows():
        cursor.execute("""
            INSERT INTO genes (gene, ensembl_id, uniprot_id, description, variant_count)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(gene) DO UPDATE SET
                ensembl_id=excluded.ensembl_id,
                uniprot_id=excluded.uniprot_id,
                description=excluded.description,
                variant_count=excluded.variant_count
        """, (
            row['Gene_Symbol'],
            row['Ensembl_Gene_ID'],
            row['UniProtKB'],
            row['Description'],
            row['Variant_Count']
        ))
    
    conn.commit()
    conn.close()
    print(f"Saved {len(genes_df)} genes")

def analyze_condition(condition="autism", max_results=50):
    """Main function to analyze a condition"""
    print(f"Analyzing: {condition}")
    
    variants = search_variants_by_condition(condition, max_results)
    if variants.empty:
        return None, None
    
    genes = create_gene_data(variants)
    
    setup_database()
    save_variants(variants, condition)
    save_genes(genes)
    
    print(f"Done: {len(variants)} variants, {len(genes)} genes")
    return variants, genes

# -------------------- DATABASE SETUP --------------------
def setup_database():
    """Initialize the SQLite database"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id TEXT,
            gene TEXT,
            clinical_significance TEXT,
            diseases TEXT,
            last_updated TEXT,
            searched_condition TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS genes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT UNIQUE,
            ensembl_id TEXT,
            uniprot_id TEXT,
            description TEXT,
            variant_count INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS analyses (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            condition TEXT,
            variant_count INTEGER,
            gene_count INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    # Insertar datos de ejemplo si la tabla está vacía
    cursor.execute("SELECT COUNT(*) FROM variants")
    if cursor.fetchone()[0] == 0:
        variants_data = [
            ("183085", "PTEN", "Pathogenic", "Autism spectrum disorder", "2023-10-15", "autism"),
            ("156468", "SHANK3", "Pathogenic", "Autism spectrum disorder", "2023-09-20", "autism"),
            ("142357", "MECP2", "Likely pathogenic", "Autism spectrum disorder", "2023-11-05", "autism"),
            ("192837", "TP53", "Pathogenic", "Breast cancer, Li-Fraumeni syndrome", "2023-12-01", "cancer"),
            ("174839", "BRCA1", "Pathogenic", "Breast cancer, Ovarian cancer", "2023-11-20", "cancer"),
            ("163849", "APOE", "Risk factor", "Alzheimer's disease", "2023-10-30", "Alzheimer's"),
            ("152739", "INS", "Pathogenic", "Diabetes mellitus", "2023-09-15", "diabetes"),
            ("142938", "CFTR", "Pathogenic", "Cystic fibrosis", "2023-11-10", "cystic fibrosis")
        ]
        
        cursor.executemany("""
            INSERT INTO variants (variant_id, gene, clinical_significance, diseases, last_updated, searched_condition)
            VALUES (?, ?, ?, ?, ?, ?)
        """, variants_data)
        
        genes_data = [
            ("PTEN", "ENSG00000171862", "P60484", "Phosphatase and tensin homolog", 3),
            ("SHANK3", "ENSG00000277372", "Q9BYB0", "SH3 and multiple ankyrin repeat domains protein 3", 2),
            ("MECP2", "ENSG00000169057", "P51608", "Methyl-CpG-binding protein 2", 1),
            ("TP53", "ENSG00000141510", "P04637", "Cellular tumor antigen p53", 5),
            ("BRCA1", "ENSG00000012048", "P38398", "Breast cancer type 1 susceptibility protein", 4),
            ("APOE", "ENSG00000130203", "P02649", "Apolipoprotein E", 2),
            ("INS", "ENSG00000254647", "P01308", "Insulin", 1),
            ("CFTR", "ENSG00000001626", "P13569", "Cystic fibrosis transmembrane conductance regulator", 1)
        ]
        
        cursor.executemany("""
            INSERT OR REPLACE INTO genes (gene, ensembl_id, uniprot_id, description, variant_count)
            VALUES (?, ?, ?, ?, ?)
        """, genes_data)
        
        print("✅ Datos de ejemplo insertados en la base de datos")
    
    conn.commit()
    conn.close()

# -------------------- RUTAS PRINCIPALES --------------------
@app.route("/")
def home():
    """Página principal"""
    setup_database()
    
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Obtener estadísticas
    cursor.execute("SELECT COUNT(*) FROM variants")
    total_variants = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genes")
    total_genes = cursor.fetchone()[0]
    
    # Obtener condiciones únicas
    cursor.execute("SELECT DISTINCT searched_condition FROM variants ORDER BY searched_condition")
    conditions = [row[0] for row in cursor.fetchall()]
    
    # Obtener análisis recientes
    cursor.execute("""
        SELECT searched_condition, COUNT(*) as variant_count, 
               COUNT(DISTINCT gene) as gene_count
        FROM variants 
        GROUP BY searched_condition 
        ORDER BY MAX(id) DESC 
        LIMIT 5
    """)
    
    recent_analyses = []
    for row in cursor.fetchall():
        recent_analyses.append({
            'condition': row[0],
            'variant_count': row[1],
            'gene_count': row[2],
            'created_at': '2024-01-01'
        })
    
    conn.close()
    
    return render_template(
        'home.html',
        total_variants=total_variants,
        total_genes=total_genes,
        conditions=conditions,
        recent_analyses=recent_analyses,
        condition_chart=None
    )

@app.route("/search", methods=['POST'])
def search():
    """Procesar búsqueda"""
    condition = request.form.get('condition', '').strip()
    max_results = int(request.form.get('max_results', 50))
    
    if condition:
        # Usar la función local analyze_condition
        variants_df, genes_df = analyze_condition(condition, max_results)
        
        # VERIFICACIÓN CORREGIDA
        if variants_df is not None and not variants_df.empty:
            # Calcular gene_count de forma segura
            gene_count = 0
            if genes_df is not None:
                gene_count = len(genes_df)
            
            # Guardar registro del análisis
            conn = sqlite3.connect(DB_FILE)
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO analyses (condition, variant_count, gene_count)
                VALUES (?, ?, ?)
            """, (condition, len(variants_df), gene_count))
            conn.commit()
            conn.close()
            
            return redirect(url_for('condition_results', condition=condition))
    
    return redirect(url_for('home'))

@app.route("/condition/<condition>")
def condition_results(condition):
    """Mostrar resultados para una condición"""
    page = request.args.get(get_page_parameter(), type=int, default=1)
    per_page = 10
    offset = (page - 1) * per_page
    
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Obtener total
    cursor.execute("SELECT COUNT(*) FROM variants WHERE searched_condition = ?", (condition,))
    total_variants = cursor.fetchone()[0]
    
    # Obtener variantes paginadas
    cursor.execute("""
        SELECT variant_id, gene, clinical_significance, diseases, last_updated
        FROM variants 
        WHERE searched_condition = ? 
        ORDER BY variant_id
        LIMIT ? OFFSET ?
    """, (condition, per_page, offset))
    
    variants_data = cursor.fetchall()
    variants = []
    for row in variants_data:
        variants.append({
            'variant_id': row[0],
            'gene': row[1],
            'clinical_significance': row[2],
            'diseases': row[3],
            'last_updated': row[4]
        })
    
    # Obtener genes para esta condición
    cursor.execute("""
        SELECT DISTINCT v.gene, g.ensembl_id, g.uniprot_id, g.description,
               COUNT(v.id) as variant_count
        FROM variants v
        LEFT JOIN genes g ON v.gene = g.gene
        WHERE v.searched_condition = ? AND v.gene IS NOT NULL
        GROUP BY v.gene
        ORDER BY variant_count DESC
    """, (condition,))
    
    genes_data = cursor.fetchall()
    genes = []
    for row in genes_data:
        genes.append({
            'gene': row[0],
            'ensembl_id': row[1],
            'uniprot_id': row[2],
            'description': row[3],
            'variant_count': row[4]
        })
    
    conn.close()
    
    pagination = Pagination(
        page=page,
        per_page=per_page,
        total=total_variants,
        record_name='variants',
        bs_version=4,
        css_framework='bootstrap4'
    )
    
    return render_template(
        'condition_results.html',
        condition=condition,
        variants=variants,
        genes=genes,
        total_variants=total_variants,
        pagination=pagination
    )

@app.route("/gene/<gene_name>")
def gene_details(gene_name):
    """Mostrar detalles de un gen"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Obtener información del gen
    cursor.execute("SELECT * FROM genes WHERE gene = ?", (gene_name,))
    gene_row = cursor.fetchone()
    
    if not gene_row:
        conn.close()
        return "Gene not found", 404
    
    gene_info = {
        'gene': gene_row[1],
        'ensembl_id': gene_row[2],
        'uniprot_id': gene_row[3],
        'description': gene_row[4],
        'variant_count': gene_row[5]
    }
    
    # Obtener variantes para este gen
    cursor.execute("""
        SELECT variant_id, clinical_significance, diseases, searched_condition, last_updated
        FROM variants 
        WHERE gene = ? 
        ORDER BY id DESC
    """, (gene_name,))
    
    variants_data = cursor.fetchall()
    variants = []
    for row in variants_data:
        variants.append({
            'variant_id': row[0],
            'clinical_significance': row[1],
            'diseases': row[2],
            'searched_condition': row[3],
            'last_updated': row[4]
        })
    
    conn.close()
    
    return render_template(
        'gene_details.html',
        gene=gene_info,
        variants=variants,
        total_variants=len(variants)
    )

@app.route('/clinvar/genes')
def clinvar_genes():
    condition = request.args.get('condition', 'autism')
    variation_type = request.args.get('variation_type', None)
    
    # Preparar filtros
    filters = {}
    if variation_type:
        filters['variation_type'] = variation_type
    
    # Consultar ClinVar con filtros
    data = query_clinvar(condition, filters)
    
    # Procesar datos
    genes = data.get('associated_genes', [])
    total_variants = data.get('total_variants', 0)
    unique_genes = data.get('unique_genes', 0)
    
    return render_template(
        'genes_list.html',
        condition=condition,
        genes=genes,
        total_variants=total_variants,
        unique_genes=unique_genes,
        filters=filters,
        current_variation_type=variation_type
    )

@app.route("/variant/<variant_id>")
def variant_details(variant_id):
    """Mostrar detalles de una variante"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT v.*, g.ensembl_id, g.uniprot_id, g.description 
        FROM variants v 
        LEFT JOIN genes g ON v.gene = g.gene 
        WHERE v.variant_id = ?
    """, (variant_id,))
    
    variant_row = cursor.fetchone()
    
    if not variant_row:
        conn.close()
        return "Variant not found", 404
    
    variant_info = {
        'variant_id': variant_row[1],
        'gene': variant_row[2],
        'clinical_significance': variant_row[3],
        'diseases': variant_row[4],
        'last_updated': variant_row[5],
        'searched_condition': variant_row[6],
        'ensembl_id': variant_row[8],
        'uniprot_id': variant_row[9],
        'description': variant_row[10]
    }
    
    cursor.execute("""
        SELECT variant_id, clinical_significance, searched_condition, diseases
        FROM variants 
        WHERE gene = ? AND variant_id != ?
        ORDER BY id DESC 
        LIMIT 10
    """, (variant_info['gene'], variant_id))
    
    similar_data = cursor.fetchall()
    similar_variants = []
    for row in similar_data:
        similar_variants.append({
            'variant_id': row[0],
            'clinical_significance': row[1],
            'searched_condition': row[2],
            'diseases': row[3]
        })
    
    conn.close()
    
    return render_template(
        'variant_details.html',
        variant=variant_info,
        similar_variants=similar_variants
    )

@app.route("/variants")
def variants_list():
    """Listar todas las variantes"""
    page = request.args.get(get_page_parameter(), type=int, default=1)
    per_page = 20
    offset = (page - 1) * per_page
    
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Obtener total
    cursor.execute("SELECT COUNT(*) FROM variants")
    total_variants = cursor.fetchone()[0]
    
    # Obtener variantes paginadas
    cursor.execute("""
        SELECT variant_id, gene, clinical_significance, diseases, searched_condition, last_updated
        FROM variants 
        ORDER BY id DESC
        LIMIT ? OFFSET ?
    """, (per_page, offset))
    
    variants_data = cursor.fetchall()
    variants = []
    for row in variants_data:
        variants.append({
            'variant_id': row[0],
            'gene': row[1],
            'clinical_significance': row[2],
            'diseases': row[3],
            'searched_condition': row[4],
            'last_updated': row[5]
        })
    
    conn.close()
    
    pagination = Pagination(
        page=page,
        per_page=per_page,
        total=total_variants,
        record_name='variants',
        bs_version=4,
        css_framework='bootstrap4'
    )
    
    return render_template(
        'variants_list.html',
        variants=variants,
        total_variants=total_variants,
        pagination=pagination
    )

@app.route("/genes")
def genes_list():
    """Listar todos los genes"""
    page = request.args.get(get_page_parameter(), type=int, default=1)
    per_page = 20
    offset = (page - 1) * per_page
    
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Obtener total
    cursor.execute("SELECT COUNT(*) FROM genes")
    total_genes = cursor.fetchone()[0]
    
    # Obtener genes paginados
    cursor.execute("""
        SELECT g.gene, g.ensembl_id, g.uniprot_id, g.description, g.variant_count,
               COUNT(v.id) as total_variants
        FROM genes g
        LEFT JOIN variants v ON g.gene = v.gene
        GROUP BY g.gene
        ORDER BY total_variants DESC
        LIMIT ? OFFSET ?
    """, (per_page, offset))
    
    genes_data = cursor.fetchall()
    genes = []
    for row in genes_data:
        genes.append({
            'gene': row[0],
            'ensembl_id': row[1],
            'uniprot_id': row[2],
            'description': row[3],
            'variant_count': row[4],
            'total_variants': row[5]
        })
    
    conn.close()
    
    pagination = Pagination(
        page=page,
        per_page=per_page,
        total=total_genes,
        record_name='genes',
        bs_version=4,
        css_framework='bootstrap4'
    )
    
    return render_template(
        'genes_list.html',
        genes=genes,
        total_genes=total_genes,
        pagination=pagination
    )

@app.route("/statistics")
def statistics():
    """Mostrar estadísticas"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Estadísticas básicas
    cursor.execute("SELECT COUNT(*) FROM variants")
    total_variants = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genes")
    total_genes = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT searched_condition) FROM variants")
    total_conditions = cursor.fetchone()[0]
    
    # Top condiciones
    cursor.execute("""
        SELECT searched_condition, COUNT(*) as count 
        FROM variants 
        GROUP BY searched_condition 
        ORDER BY count DESC 
        LIMIT 15
    """)
    
    top_conditions_data = cursor.fetchall()
    top_conditions = []
    for row in top_conditions_data:
        top_conditions.append({
            'searched_condition': row[0],
            'count': row[1]
        })
    
    # Top genes
    cursor.execute("""
        SELECT gene, COUNT(*) as count 
        FROM variants 
        WHERE gene IS NOT NULL AND gene != ''
        GROUP BY gene 
        ORDER BY count DESC 
        LIMIT 15
    """)
    
    top_genes_data = cursor.fetchall()
    top_genes = []
    for row in top_genes_data:
        top_genes.append({
            'gene': row[0],
            'count': row[1]
        })
    
    # Distribución de significancia clínica
    cursor.execute("""
        SELECT clinical_significance, COUNT(*) as count
        FROM variants 
        WHERE clinical_significance IS NOT NULL AND clinical_significance != ''
        GROUP BY clinical_significance 
        ORDER BY count DESC
        LIMIT 10
    """)
    
    significance_data = cursor.fetchall()
    significance_dist = []
    for row in significance_data:
        significance_dist.append({
            'clinical_significance': row[0],
            'count': row[1]
        })
    
    conn.close()
    
    return render_template(
        'statistics.html',
        total_variants=total_variants,
        total_genes=total_genes,
        total_conditions=total_conditions,
        top_conditions=top_conditions,
        top_genes=top_genes,
        significance_dist=significance_dist,
        charts={}  # Sin gráficos por ahora
    )

@app.route("/api/genes_list")
def api_genes_list():
    """API para autocomplete de genes"""
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT gene FROM genes WHERE gene IS NOT NULL ORDER BY gene")
    genes = [row[0] for row in cursor.fetchall()]
    conn.close()
    
    return jsonify({"genes": genes})

# -------------------- ERROR HANDLERS --------------------
@app.errorhandler(404)
def page_not_found(e):
    return render_template('error.html', message="Page not found"), 404

@app.errorhandler(500)
def internal_error(e):
    return render_template('error.html', message="Internal server error"), 500

if __name__ == "__main__":
    setup_database()
    app.run(debug=True, port=5001, host='0.0.0.0')
