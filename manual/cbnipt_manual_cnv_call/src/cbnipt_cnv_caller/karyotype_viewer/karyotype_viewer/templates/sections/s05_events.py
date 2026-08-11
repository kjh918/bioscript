"""
S05 — Chromosomal Events
이상 염색체 이벤트 카드 — ISCN 크게, 타입/chr/CN 표시
JS: renderEvCards()
"""

SECTION_EVENTS = """
<!-- ═══════════════════════════════════════════════════════ S05 EVENTS -->
<div class="card" style="margin-bottom:1rem;">
  <div class="ch">
    <span class="ct">Chromosomal Events</span>
    <span style="font-size:10px;color:var(--text-muted);">검출된 염색체 이상</span>
  </div>
  <div class="cb" style="padding:12px 14px;">
    <div class="ev-cards" id="ev-cards"></div>
  </div>
</div>
"""