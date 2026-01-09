#!/usr/bin/env python3
"""
Generate missing embeddings for Mirror7 Knowledge Base.
Uses sentence-transformers with all-MiniLM-L6-v2 (384 dimensions).
"""

import json
import urllib.request
import urllib.error
import base64
import sys

try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    print("❌ sentence-transformers non installato!")
    print("   Esegui: pip install sentence-transformers")
    sys.exit(1)

# Config
SURREAL_URL = "http://localhost:8000/sql"
SURREAL_NS = "research"
SURREAL_DB = "knowledge"
SURREAL_AUTH = base64.b64encode(b"root:root").decode()
BATCH_SIZE = 50

# Load model
print("🔄 Loading embedding model...")
model = SentenceTransformer('all-MiniLM-L6-v2')
print("✅ Model loaded (384 dimensions)")

def query_surreal(sql: str):
    """Execute SQL query against SurrealDB."""
    data = sql.encode('utf-8')
    req = urllib.request.Request(
        SURREAL_URL,
        data=data,
        headers={
            'Accept': 'application/json',
            'Content-Type': 'text/plain',
            'Authorization': f'Basic {SURREAL_AUTH}',
            'surreal-ns': SURREAL_NS,
            'surreal-db': SURREAL_DB,
        },
        method='POST'
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(f"❌ Query error: {e}")
        return []

def get_stats():
    """Get embedding statistics."""
    result = query_surreal("""
        SELECT count() as total FROM knowledge GROUP ALL;
    """)
    total = result[0]['result'][0]['total'] if result else 0
    
    result = query_surreal("""
        SELECT count() as with_emb FROM knowledge 
        WHERE embedding != NONE AND array::len(embedding) > 0 GROUP ALL;
    """)
    with_emb = result[0]['result'][0]['with_emb'] if result and result[0].get('result') else 0
    
    return total, with_emb

def get_missing(limit: int = BATCH_SIZE):
    """Get records without embeddings."""
    result = query_surreal(f"""
        SELECT id, title, content FROM knowledge 
        WHERE embedding = NONE OR array::len(embedding) = 0 
        LIMIT {limit};
    """)
    if result and result[0].get('result'):
        return result[0]['result']
    return []

def generate_embedding(text: str):
    """Generate 384-dim embedding."""
    return model.encode(text, normalize_embeddings=True).tolist()

def update_embedding(record_id: str, embedding: list) -> bool:
    """Update record with embedding."""
    # Format as SurrealQL array
    emb_str = "[" + ",".join(f"{v}f" for v in embedding) + "]"
    sql = f"UPDATE {record_id} SET embedding = {emb_str};"
    result = query_surreal(sql)
    return bool(result)

def main():
    print("\n" + "=" * 60)
    print("🧠 MIRROR7 EMBEDDING GENERATOR")
    print("=" * 60)
    
    # Initial stats
    total, with_emb = get_stats()
    missing = total - with_emb
    print(f"\n📊 INITIAL STATE:")
    print(f"   Total: {total}")
    print(f"   With embeddings: {with_emb}")
    print(f"   Missing: {missing}")
    
    if missing == 0:
        print("\n✅ All records already have embeddings!")
        return
    
    print(f"\n🚀 Generating embeddings for {missing} records...")
    print(f"   Batch size: {BATCH_SIZE}")
    
    generated = 0
    errors = 0
    batch_num = 0
    
    while True:
        records = get_missing(BATCH_SIZE)
        if not records:
            break
        
        batch_num += 1
        print(f"\n📦 Batch {batch_num} ({len(records)} records)")
        
        for rec in records:
            try:
                # Combine title + content for richer embedding
                title = rec.get('title', '') or ''
                content = rec.get('content', '') or ''
                text = f"{title}\n{content}"[:8000]  # Limit length
                
                if not text.strip():
                    errors += 1
                    continue
                
                # Generate and update
                embedding = generate_embedding(text)
                if update_embedding(str(rec['id']), embedding):
                    generated += 1
                    print(f"   ✓ {title[:50]}...")
                else:
                    errors += 1
                    
            except Exception as e:
                errors += 1
                print(f"   ✗ Error: {e}")
    
    # Final stats
    total, with_emb = get_stats()
    print(f"\n📊 FINAL STATE:")
    print(f"   Total: {total}")
    print(f"   With embeddings: {with_emb}")
    print(f"   Coverage: {with_emb/total*100:.1f}%")
    
    print(f"\n✅ DONE: {generated} embeddings generated, {errors} errors")

if __name__ == "__main__":
    main()
