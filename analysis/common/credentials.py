"""Credential helpers.

No tokens are stored in this repository. Set them in the environment:

    export CAVE_TOKEN="..."
    export CATMAID_API_TOKEN="..."

Notebooks in this repo that previously carried an inline token now contain the literal
placeholder ``REPLACE_WITH_YOUR_CAVE_TOKEN`` / ``REPLACE_WITH_YOUR_CATMAID_API_TOKEN``.
Either substitute your own token there, or replace that cell with a call into this module.
"""
import os

CAVE_DATASTACK = "wclee_mouse_spinalcord_cltmr"
SEGMENTATION_SOURCE = (
    "graphene://https://cave.fanc-fly.com/segmentation/table/wclee_mouse_spinalcord_cltmr"
)
CATMAID_SERVER = "https://radagast.hms.harvard.edu/dorsalhorn/"


def cave_token():
    """CAVE token from $CAVE_TOKEN, else from the standard caveclient credentials file."""
    token = os.environ.get("CAVE_TOKEN")
    if token:
        return token
    import caveclient

    token = caveclient.auth.AuthClient().token
    if not token:
        raise RuntimeError(
            "No CAVE token. Set $CAVE_TOKEN, or run:\n"
            "  import caveclient; "
            "caveclient.auth.AuthClient().save_token(token='your-token')"
        )
    return token


def cave_client():
    """A CAVEclient for the dSC1 datastack."""
    import caveclient

    return caveclient.CAVEclient(CAVE_DATASTACK, auth_token=cave_token())


def cloud_volume():
    """A CloudVolume handle on the dSC1 segmentation."""
    import cloudvolume as cv

    return cv.CloudVolume(
        SEGMENTATION_SOURCE, secrets=cave_token(), use_https=True, progress=False
    )


def catmaid_token():
    token = os.environ.get("CATMAID_API_TOKEN")
    if not token:
        raise RuntimeError("Set $CATMAID_API_TOKEN to your CATMAID API token.")
    return token


def catmaid_instance(project_id, http_user=None, http_password=None):
    """A pymaid CatmaidInstance for one of the dSC_APEX projects."""
    import pymaid

    return pymaid.CatmaidInstance(
        server=CATMAID_SERVER,
        api_token=catmaid_token(),
        http_user=http_user or os.environ.get("CATMAID_HTTP_USER"),
        http_password=http_password or os.environ.get("CATMAID_HTTP_PASSWORD"),
        project_id=project_id,
    )
