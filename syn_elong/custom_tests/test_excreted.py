from memote.utils import annotate, wrapper
from syn_elong import model as syn, expected_metab
from concerto.testing.secretion import get_excreted_metabolites


@annotate(
    title="Models ability to excrete metabolites",
    format_type="number"
)
def test_excreted():
    """
    Tests models ability to excrete metabolites

    This test addresses the models ability to excrete metabolites.
    Models are optimized to produce the given exchange reaction and
    simulated with natural medium.
    """
    ann = test_excreted.annotation
    tp, fn = get_excreted_metabolites(syn, expected_metab)
    ann["data"] = {'TP': sorted(tp), 'FN': sorted(fn)}
    ann["metric"] = len(ann["data"]['FN'])
    ann["message"] = wrapper.fill(
        f"""Model is able to excrete {len(ann['data']['TP'])} and unable to 
        excrete {len(ann['data']['FN'])}"""
        )
    assert len(ann['data']['FN']) == 0, ann["message"]
