from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from gbdraw.features.objects import FeatureLocationPart, FeatureObject


def make_origin_spanning_seq_feature() -> SeqFeature:
    return SeqFeature(
        CompoundLocation(
            [
                SimpleLocation(0, 576, strand=-1),
                SimpleLocation(16023, 16569, strand=-1),
            ],
            operator="join",
        ),
        type="D-loop",
        qualifiers={},
    )


def make_origin_spanning_feature_object(record_id: str = "rec1") -> FeatureObject:
    return FeatureObject(
        feature_id="feature_000000099",
        location=[
            FeatureLocationPart("block", "001", "negative", 1, 576, False),
            FeatureLocationPart("block", "002", "negative", 16023, 16569, True),
        ],
        is_directional=False,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=[
            SimpleLocation(0, 576, strand=-1),
            SimpleLocation(16023, 16569, strand=-1),
        ],
        type="D-loop",
        qualifiers={},
        record_id=record_id,
    )
