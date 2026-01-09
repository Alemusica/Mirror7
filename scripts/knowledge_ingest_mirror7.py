#!/usr/bin/env python3
"""
Mirror7 Knowledge Ingestion System

Extended version of knowledge_ingest.py with:
- Binaural/neuroscience specific queries
- PubMed support
- Semantic Scholar API
- Non-conventional source queries
- Golden ratio / Fibonacci research

Usage:
    python knowledge_ingest_mirror7.py --mode binaural    # Core binaural research
    python knowledge_ingest_mirror7.py --mode neuro       # Neuroscience
    python knowledge_ingest_mirror7.py --mode phi         # Golden ratio research
    python knowledge_ingest_mirror7.py --mode alternative # Non-conventional
    python knowledge_ingest_mirror7.py --mode all         # Everything
    python knowledge_ingest_mirror7.py --stats            # Show stats
"""

import json
import hashlib
import argparse
import urllib.request
import urllib.error
import urllib.parse
import time
import ssl
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import List, Dict, Set, Optional

# SSL context for macOS
ssl_context = ssl.create_default_context()
ssl_context.check_hostname = False
ssl_context.verify_mode = ssl.CERT_NONE

# ============================================================================
# CONFIGURATION
# ============================================================================

SURREALDB_URL = "http://localhost:8000/sql"
SURREALDB_NS = "research"
SURREALDB_DB = "knowledge"

ARXIV_API = "http://export.arxiv.org/api/query"
PUBMED_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SEMANTIC_SCHOLAR_API = "https://api.semanticscholar.org/graph/v1"

# ============================================================================
# QUERY SETS FOR MIRROR7
# ============================================================================

QUERY_SETS = {
    "binaural": {
        "description": "Core binaural and spatial audio research",
        "arxiv": [
            "binaural perception interaural time difference",
            "ITD ILD sound localization neural",
            "head related transfer function HRTF",
            "spatial audio 3D sound",
            "binaural rendering virtual acoustics",
            "auditory spatial perception",
        ],
        "pubmed": [
            "binaural hearing localization",
            "interaural time difference perception",
            "auditory spatial processing brain",
        ],
        "tags": ["binaural", "spatial-audio", "itd", "ild", "hrtf"],
    },
    
    "neuro": {
        "description": "Neuroscience of auditory processing and hemispheric communication",
        "arxiv": [
            "auditory cortex hemispheric lateralization",
            "corpus callosum interhemispheric transfer",
            "binaural beat entrainment EEG",
            "alpha wave synchronization auditory",
            "neural oscillation frequency entrainment",
            "dichotic listening hemispheric",
            # NEW: User requested queries
            "alpha wave binaural entrainment EEG",
            "corpus callosum auditory transfer",
            "dichotic listening hemispheric lateralization",
            "40 Hz gamma entrainment Alzheimer",
            "theta wave meditation auditory",
            "brain hemisphere synchronization sound",
        ],
        "pubmed": [
            "binaural beats brain hemispheric synchronization",
            "corpus callosum auditory processing",
            "interhemispheric coherence EEG",
            "auditory cortex lateralization fMRI",
            "alpha wave entrainment meditation",
            # NEW: Clinical/medical queries
            "gamma oscillation 40Hz cognitive",
            "dichotic listening schizophrenia",
            "binaural beats anxiety treatment",
            "auditory steady state response ASSR",
        ],
        "tags": ["neuroscience", "eeg", "hemispheric", "entrainment", "cortex", "gamma", "alpha"],
    },
    
    "phi": {
        "description": "Golden ratio, Fibonacci, and mathematical beauty in sound",
        "arxiv": [
            "golden ratio acoustics vibration",
            "Fibonacci sequence music harmony",
            "phi proportion natural patterns",
            "golden angle phyllotaxis spiral",
            "irrational numbers quasi-periodic",
            "Penrose tiling aperiodic",
        ],
        "pubmed": [
            "golden ratio aesthetic perception",
            "Fibonacci biological patterns",
        ],
        "tags": ["golden-ratio", "phi", "fibonacci", "mathematics", "patterns"],
    },
    
    "alternative": {
        "description": "Non-conventional frequency and consciousness research",
        "arxiv": [
            "cymatics pattern formation vibration",
            "Schumann resonance biological effects",
            "infrasound physiological response",
            "heart rate variability coherence",
            "low frequency sound health",
            "vibroacoustic therapy",
        ],
        "pubmed": [
            "frequency therapy systematic review",
            "sound healing clinical trial",
            "heart rate variability biofeedback",
            "Schumann resonance human health",
            "infrasound effects human body",
            "music therapy neurological",
        ],
        "tags": ["frequency-therapy", "cymatics", "hrv", "alternative", "vibroacoustic"],
    },
    
    "russian": {
        "description": "Soviet/Russian bioacoustics and unconventional research",
        "arxiv": [
            "Russian bioacoustics research",
            "Soviet infrasound military",
            "wave genetics DNA",
            "torsion field physics",
            "biophoton cell communication",
            "millimeter wave therapy",
        ],
        "pubmed": [
            "Soviet biofeedback research",
            "Russian psychophysiology",
            "millimeter wave biological effects",
            "biophoton emission cells",
        ],
        "tags": ["russian", "soviet", "bioacoustics", "unconventional"],
    },
    
    "corpus_callosum": {
        "description": "Corpus callosum stimulation and hemispheric integration",
        "arxiv": [
            "corpus callosum auditory stimulation",
            "transcallosal transfer auditory",
            "interhemispheric communication sound",
            "whole brain synchronization",
        ],
        "pubmed": [
            "corpus callosum stimulation therapy",
            "transcallosal auditory processing",
            "hemispheric integration binaural",
            "interhemispheric coherence EEG sound",
            "agenesis corpus callosum auditory",
        ],
        "tags": ["corpus-callosum", "hemispheric", "integration", "transcallosal"],
    },
    
    "gamma_therapy": {
        "description": "40 Hz gamma stimulation and cognitive enhancement",
        "arxiv": [
            "40 Hz gamma entrainment cognition",
            "gamma oscillation stimulation brain",
            "sensory stimulation gamma Alzheimer",
        ],
        "pubmed": [
            "40 Hz gamma stimulation Alzheimer",
            "gamma entrainment cognitive enhancement",
            "auditory steady state response 40Hz",
            "gamma oscillation memory",
            "light sound therapy gamma",
        ],
        "tags": ["gamma", "40hz", "entrainment", "cognitive", "therapy"],
    },
    
    "heart_brain": {
        "description": "Heart-brain coherence and HRV biofeedback",
        "arxiv": [
            "heart brain coherence",
            "HRV biofeedback neural",
        ],
        "pubmed": [
            "HeartMath coherence research",
            "heart rate variability biofeedback",
            "vagal tone meditation",
            "resonance frequency breathing",
            "cardiac coherence training",
        ],
        "tags": ["hrv", "heart", "coherence", "biofeedback", "vagal"],
    },
    
    "implementation": {
        "description": "Practical DSP implementation guides and code",
        "arxiv": [
            "binaural synthesis algorithm",
            "granular synthesis real-time",
            "HRTF convolution efficient",
            "spatial audio DSP",
        ],
        "pubmed": [],  # Not relevant for implementation
        "stackexchange": [
            # DSP StackExchange
            ("dsp", "binaural panning algorithm"),
            ("dsp", "ITD ILD implementation"),
            ("dsp", "HRTF convolution"),
            ("dsp", "granular synthesis"),
            ("dsp", "pink noise generation"),
            ("dsp", "golden ratio filter"),
            # Audio StackExchange
            ("sound", "binaural audio creation"),
            ("sound", "spatial audio mixing"),
            # Signal Processing
            ("dsp", "delay line interpolation"),
            ("dsp", "allpass filter reverb"),
            ("dsp", "comb filter design"),
        ],
        "tags": ["implementation", "dsp", "algorithm", "code"],
    },
}

# ============================================================================
# SURREALDB HELPERS (copied from original)
# ============================================================================

import base64
SURREALDB_AUTH = base64.b64encode(b"root:root").decode()

def surreal_query(sql: str, silent: bool = False) -> dict:
    """Execute SQL query against SurrealDB."""
    data = sql.encode('utf-8')
    req = urllib.request.Request(
        SURREALDB_URL,
        data=data,
        headers={
            'Accept': 'application/json',
            'Content-Type': 'text/plain',
            'Authorization': f'Basic {SURREALDB_AUTH}',
            'surreal-ns': SURREALDB_NS,
            'surreal-db': SURREALDB_DB,
        },
        method='POST'
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        if not silent:
            print(f"SurrealDB Error: {e}")
        return {}

def escape_surreal(s: str) -> str:
    """Escape string for SurrealDB."""
    if not s:
        return ""
    s = s.replace('\\', '\\\\')
    s = s.replace('"', '\\"')
    s = s.replace('\n', '\\n')
    s = s.replace('\r', '\\r')
    s = s.replace('\t', '\\t')
    return s

def get_existing_ids(source: str) -> Set[str]:
    """Get already ingested IDs to avoid duplicates."""
    sql = f'SELECT source_id FROM knowledge WHERE source = "{source}";'
    result = surreal_query(sql, silent=True)
    if result and result[0].get('result'):
        return {r['source_id'] for r in result[0]['result']}
    return set()

# ============================================================================
# ARXIV INGEST (improved)
# ============================================================================

def fetch_arxiv(query: str, start: int = 0, max_results: int = 50) -> List[Dict]:
    """Fetch papers from ArXiv API with rate limiting."""
    params = {
        'search_query': f'all:{query}',
        'start': start,
        'max_results': max_results,
        'sortBy': 'relevance',
        'sortOrder': 'descending',
    }
    
    url = f"{ARXIV_API}?{urllib.parse.urlencode(params)}"
    
    # ArXiv rate limit: 1 request per 3 seconds
    time.sleep(3)
    
    try:
        with urllib.request.urlopen(url, timeout=30, context=ssl_context) as resp:
            xml_data = resp.read().decode()
            
        root = ET.fromstring(xml_data)
        ns = {'atom': 'http://www.w3.org/2005/Atom'}
        
        papers = []
        for entry in root.findall('atom:entry', ns):
            paper = {
                'id': entry.find('atom:id', ns).text.split('/')[-1] if entry.find('atom:id', ns) is not None else '',
                'title': entry.find('atom:title', ns).text.strip().replace('\n', ' ') if entry.find('atom:title', ns) is not None else '',
                'summary': entry.find('atom:summary', ns).text.strip() if entry.find('atom:summary', ns) is not None else '',
                'published': entry.find('atom:published', ns).text if entry.find('atom:published', ns) is not None else '',
                'authors': [a.find('atom:name', ns).text for a in entry.findall('atom:author', ns) if a.find('atom:name', ns) is not None],
                'categories': [c.get('term') for c in entry.findall('atom:category', ns)],
                'link': entry.find('atom:id', ns).text if entry.find('atom:id', ns) is not None else '',
            }
            papers.append(paper)
        
        return papers
    except Exception as e:
        print(f"  ArXiv API error: {e}")
        return []

def ingest_arxiv_query(query: str, tags: List[str], max_papers: int = 20) -> int:
    """Ingest ArXiv papers for a specific query."""
    print(f"  📄 ArXiv: '{query}'")
    
    existing_ids = get_existing_ids('arxiv')
    papers = fetch_arxiv(query, max_results=max_papers + 10)
    ingested = 0
    
    for paper in papers:
        if ingested >= max_papers:
            break
        
        arxiv_id_raw = paper['id']
        if arxiv_id_raw in existing_ids:
            continue
        
        arxiv_id = arxiv_id_raw.replace('.', '_').replace('/', '_')
        title = escape_surreal(paper['title'])
        content = escape_surreal(paper['summary'][:10000])
        url = paper['link']
        
        published = paper.get('published', '')
        created_iso = None
        if published:
            try:
                dt = datetime.fromisoformat(published.replace('Z', '+00:00'))
                created_iso = dt.strftime('%Y-%m-%dT%H:%M:%SZ')
            except:
                pass
        
        metadata = {
            'authors': paper.get('authors', []),
            'categories': paper.get('categories', []),
            'query': query,
        }
        
        all_tags = list(set(tags + paper.get('categories', [])))
        
        sql = f'''
        UPSERT knowledge:arxiv_{arxiv_id} CONTENT {{
            source: "arxiv",
            source_id: "{escape_surreal(arxiv_id_raw)}",
            title: "{title}",
            content: "{content}",
            url: "{escape_surreal(url)}",
            created_at: {f'd"{created_iso}"' if created_iso else 'NONE'},
            metadata: {json.dumps(metadata)},
            tags: {json.dumps(all_tags)}
        }};
        '''
        result = surreal_query(sql, silent=True)
        if result:
            ingested += 1
            print(f"    ✓ {title[:50]}...")
    
    return ingested

# ============================================================================
# PUBMED INGEST
# ============================================================================

def fetch_pubmed(query: str, max_results: int = 20) -> List[Dict]:
    """Fetch papers from PubMed API."""
    # Step 1: Search for IDs
    search_params = {
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json',
        'sort': 'relevance',
    }
    
    search_url = f"{PUBMED_API}/esearch.fcgi?{urllib.parse.urlencode(search_params)}"
    
    try:
        with urllib.request.urlopen(search_url, timeout=30, context=ssl_context) as resp:
            search_result = json.loads(resp.read().decode())
        
        ids = search_result.get('esearchresult', {}).get('idlist', [])
        if not ids:
            return []
        
        # Step 2: Fetch details
        time.sleep(0.5)  # Rate limit
        
        fetch_params = {
            'db': 'pubmed',
            'id': ','.join(ids),
            'retmode': 'xml',
        }
        
        fetch_url = f"{PUBMED_API}/efetch.fcgi?{urllib.parse.urlencode(fetch_params)}"
        
        with urllib.request.urlopen(fetch_url, timeout=30, context=ssl_context) as resp:
            xml_data = resp.read().decode()
        
        root = ET.fromstring(xml_data)
        
        papers = []
        for article in root.findall('.//PubmedArticle'):
            pmid = article.find('.//PMID')
            title = article.find('.//ArticleTitle')
            abstract = article.find('.//AbstractText')
            
            paper = {
                'id': pmid.text if pmid is not None else '',
                'title': title.text if title is not None else '',
                'abstract': abstract.text if abstract is not None else '',
                'link': f"https://pubmed.ncbi.nlm.nih.gov/{pmid.text}/" if pmid is not None else '',
            }
            papers.append(paper)
        
        return papers
    except Exception as e:
        print(f"  PubMed API error: {e}")
        return []

def ingest_pubmed_query(query: str, tags: List[str], max_papers: int = 10) -> int:
    """Ingest PubMed papers for a specific query."""
    print(f"  🏥 PubMed: '{query}'")
    
    existing_ids = get_existing_ids('pubmed')
    papers = fetch_pubmed(query, max_results=max_papers + 5)
    ingested = 0
    
    for paper in papers:
        if ingested >= max_papers:
            break
        
        pmid = paper['id']
        if not pmid or pmid in existing_ids:
            continue
        
        title = escape_surreal(paper.get('title') or 'Untitled')
        abstract_text = paper.get('abstract') or ''
        content = escape_surreal(abstract_text[:10000])
        url = paper['link']
        
        metadata = {'query': query}
        
        sql = f'''
        UPSERT knowledge:pubmed_{pmid} CONTENT {{
            source: "pubmed",
            source_id: "{pmid}",
            title: "{title}",
            content: "{content}",
            url: "{escape_surreal(url)}",
            metadata: {json.dumps(metadata)},
            tags: {json.dumps(tags)}
        }};
        '''
        result = surreal_query(sql, silent=True)
        if result:
            ingested += 1
            print(f"    ✓ {title[:50]}...")
    
    return ingested

# ============================================================================
# STACKEXCHANGE INGEST
# ============================================================================

STACKEXCHANGE_API = "https://api.stackexchange.com/2.3"

def fetch_stackexchange(site: str, query: str, max_results: int = 15) -> List[Dict]:
    """Fetch Q&A from StackExchange API."""
    params = {
        'order': 'desc',
        'sort': 'relevance',
        'intitle': query,
        'site': site,
        'filter': 'withbody',
        'pagesize': min(max_results, 30),
    }
    
    url = f"{STACKEXCHANGE_API}/search/advanced?{urllib.parse.urlencode(params)}"
    
    try:
        req = urllib.request.Request(url, headers={'Accept-Encoding': 'gzip'})
        with urllib.request.urlopen(req, timeout=15, context=ssl_context) as resp:
            import gzip
            if resp.info().get('Content-Encoding') == 'gzip':
                data = gzip.decompress(resp.read())
            else:
                data = resp.read()
            result = json.loads(data.decode())
        
        questions = []
        for item in result.get('items', []):
            questions.append({
                'id': str(item.get('question_id', '')),
                'title': item.get('title', ''),
                'body': item.get('body', '')[:5000],  # Truncate
                'link': item.get('link', ''),
                'score': item.get('score', 0),
                'tags': item.get('tags', []),
            })
        return questions
    except Exception as e:
        print(f"    SE API error ({site}): {e}")
        return []

def ingest_stackexchange_query(site: str, query: str, tags: List[str], max_items: int = 10) -> int:
    """Ingest StackExchange Q&A for a specific query."""
    print(f"  💬 {site}.SE: '{query}'")
    
    source_name = f"{site}.stackexchange"
    existing_ids = get_existing_ids(source_name)
    questions = fetch_stackexchange(site, query, max_results=max_items + 5)
    ingested = 0
    
    for q in questions:
        if ingested >= max_items:
            break
        
        qid = q['id']
        if not qid or qid in existing_ids:
            continue
        
        title = escape_surreal(q.get('title') or 'Untitled')
        # Strip HTML tags from body
        import re
        body = re.sub(r'<[^>]+>', '', q.get('body', ''))
        content = escape_surreal(body[:8000])
        url = q.get('link', '')
        
        metadata = {
            'score': q.get('score', 0),
            'query': query,
            'question_tags': q.get('tags', []),
        }
        
        all_tags = list(set(tags + q.get('tags', [])))
        
        sql = f'''
        UPSERT knowledge:{site}_{qid} CONTENT {{
            source: "{source_name}",
            source_id: "{qid}",
            title: "{title}",
            content: "{content}",
            url: "{escape_surreal(url)}",
            metadata: {json.dumps(metadata)},
            tags: {json.dumps(all_tags)}
        }};
        '''
        result = surreal_query(sql, silent=True)
        if result:
            ingested += 1
            print(f"    ✓ [{q.get('score',0):+d}] {title[:45]}...")
    
    return ingested

# ============================================================================
# MAIN INGESTION LOGIC
# ============================================================================

def ingest_query_set(mode: str, max_per_query: int = 10):
    """Ingest all queries for a specific mode."""
    if mode not in QUERY_SETS:
        print(f"❌ Unknown mode: {mode}")
        print(f"   Available: {', '.join(QUERY_SETS.keys())}")
        return
    
    query_set = QUERY_SETS[mode]
    tags = query_set['tags']
    
    print(f"\n{'='*60}")
    print(f"📚 {query_set['description'].upper()}")
    print(f"{'='*60}")
    
    total_ingested = 0
    
    # ArXiv queries
    arxiv_queries = query_set.get('arxiv', [])
    if arxiv_queries:
        print(f"\n🔬 ArXiv ({len(arxiv_queries)} queries)")
        for query in arxiv_queries:
            count = ingest_arxiv_query(query, tags, max_papers=max_per_query)
            total_ingested += count
    
    # PubMed queries
    pubmed_queries = query_set.get('pubmed', [])
    if pubmed_queries:
        print(f"\n🏥 PubMed ({len(pubmed_queries)} queries)")
        for query in pubmed_queries:
            count = ingest_pubmed_query(query, tags, max_papers=max_per_query // 2)
            total_ingested += count
    
    # StackExchange queries
    se_queries = query_set.get('stackexchange', [])
    if se_queries:
        print(f"\n💬 StackExchange ({len(se_queries)} queries)")
        for site, query in se_queries:
            count = ingest_stackexchange_query(site, query, tags, max_items=max_per_query // 2)
            total_ingested += count
            time.sleep(1)  # Rate limit for SE API
    
    print(f"\n✅ Total ingested for '{mode}': {total_ingested} items")
    return total_ingested

def show_stats():
    """Show knowledge base statistics."""
    print("\n" + "=" * 60)
    print("MIRROR7 KNOWLEDGE BASE STATISTICS")
    print("=" * 60)
    
    result = surreal_query("SELECT source, count() FROM knowledge GROUP BY source;")
    if result and result[0].get('result'):
        print("\n📊 By Source:")
        for row in result[0]['result']:
            print(f"   {row.get('source', 'unknown')}: {row.get('count', 0)} items")
    
    result = surreal_query("SELECT count() FROM knowledge GROUP ALL;")
    if result and result[0].get('result'):
        total = result[0]['result'][0].get('count', 0)
        print(f"\n   TOTAL: {total} knowledge items")
    
    # By tags
    result = surreal_query("""
        SELECT tags, count() FROM knowledge 
        WHERE tags != NONE 
        GROUP BY tags 
        ORDER BY count DESC 
        LIMIT 20;
    """, silent=True)
    if result and result[0].get('result'):
        print("\n🏷️ Top Tags:")
        for row in result[0]['result'][:10]:
            tags = row.get('tags', [])
            if tags:
                print(f"   {tags}: {row.get('count', 0)}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Mirror7 Knowledge Ingestion System')
    parser.add_argument('--mode', 
                        choices=['binaural', 'neuro', 'phi', 'alternative', 'russian', 
                                 'corpus_callosum', 'gamma_therapy', 'heart_brain',
                                 'implementation', 'all'],
                        default='all', help='Query set to ingest')
    parser.add_argument('--max', type=int, default=15, help='Max papers per query')
    parser.add_argument('--stats', action='store_true', help='Show KB statistics')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("MIRROR7 KNOWLEDGE INGESTION SYSTEM")
    print("Binaural • Neuroscience • Φ Research • Alternative")
    print("=" * 60)
    
    # Test connection
    print("\n🔌 Connecting to SurrealDB...")
    result = surreal_query("INFO FOR DB;", silent=True)
    if not result:
        print("❌ Cannot connect to SurrealDB!")
        print("Start it with:")
        print("  surreal start --log warn --user root --pass root file:~/.config/surrealdb/knowledge.db")
        return
    print("✅ Connected!")
    
    if args.stats:
        show_stats()
        return
    
    # Ingest based on mode
    if args.mode == 'all':
        for mode in QUERY_SETS.keys():
            ingest_query_set(mode, max_per_query=args.max)
    else:
        ingest_query_set(args.mode, max_per_query=args.max)
    
    # Show final stats
    show_stats()
    
    print("\n📖 Available modes:")
    for mode, config in QUERY_SETS.items():
        print(f"   --mode {mode:12} → {config['description']}")

if __name__ == "__main__":
    main()
