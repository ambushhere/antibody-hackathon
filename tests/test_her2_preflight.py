from her2_preflight import annotate
from her2_candidates import VH_PREFIX,VH_SUFFIX,TRASTUZUMAB_VL


def test_wt_is_numberable_but_not_novel():
    r=annotate(VH_PREFIX+'WGGDGFYAMD'+VH_SUFFIX,TRASTUZUMAB_VL)
    assert r['numbering_pass']
    assert r['canonical_cysteines_only']
    assert not r['all_schemes_identity_below_70']
    assert r['cdrh3_by_scheme']['imgt']['sequence']=='SRWGGDGFYAMDY'
    assert r['cdrh3_by_scheme']['kabat']['sequence']=='WGGDGFYAMDY'


def test_extra_cysteine_fails_local_gate():
    r=annotate(VH_PREFIX+'CGGDGFYAMD'+VH_SUFFIX,TRASTUZUMAB_VL)
    assert not r['canonical_cysteines_only']
