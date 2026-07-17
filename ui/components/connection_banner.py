"""One banner for "the server is not answering", instead of ten identical toasts.

When the Dash server stops responding — on this box, usually because the OOM killer
took it — every in-flight callback fails independently, and dash-renderer stacks one
"Callback failed: the server did not respond" toast per failure. Ten identical red
cards look like ten separate catastrophes; they are one, repeated.

This collapses that into a single banner carrying a retry count, and it says something
actionable, because on a memory-starved machine the cause is almost never the callback
that happened to be in flight.
"""

from __future__ import annotations

from dash import Input, Output, State, dcc, html

BANNER_ID = "connection-banner"
TEXT_ID = "connection-banner-text"
PROBE_ID = "connection-probe"
STATE_ID = "connection-state"

_HIDDEN = {"display": "none"}
_VISIBLE = {
    "display": "block",
    "position": "sticky",
    "top": "0",
    "zIndex": "1000",
    "padding": "10px 14px",
    "borderRadius": "8px",
    "margin": "0 0 12px 0",
    "background": "rgba(248,81,73,0.12)",
    "border": "1px solid rgba(248,81,73,0.45)",
    "color": "#f85149",
    "fontWeight": "600",
}


DETAIL_ID = "connection-banner-detail"


def connection_banner():
    """The banner plus the heartbeat that drives it. Mount once, in the app layout."""
    return html.Div(
        [
            # 5s notices a dead server quickly and costs one tiny JSON response.
            dcc.Interval(id=PROBE_ID, interval=5000, n_intervals=0),
            dcc.Store(id=STATE_ID, data={"failures": 0, "first": None, "last": None}),
            html.Div(
                [
                    html.Span(id=TEXT_ID),
                    # The individual occurrences, deduplicated, for anyone who wants them.
                    # One banner + a disclosure beats ten identical red cards.
                    html.Details(
                        [
                            html.Summary("Ver ocorrências", className="details-summary", style={"padding": "6px 0"}),
                            html.Div(id=DETAIL_ID, className="mini-note", style={"fontFamily": "monospace"}),
                        ],
                        className="details-panel",
                        style={"background": "transparent", "border": "none", "marginTop": "6px"},
                    ),
                ],
                id=BANNER_ID,
                style=_HIDDEN,
            ),
        ]
    )


def register_clientside(app) -> None:
    """Poll /health from the browser and drive the banner.

    Both callbacks must be clientside: a server-side callback cannot report that the
    server is down, because running it needs the very server that is down.
    """
    app.clientside_callback(
        """
        function(n, state) {
            const prev  = (state && state.failures) || 0;
            const first = (state && state.first) || null;
            const now   = new Date().toTimeString().slice(0, 8);
            return fetch('/health', {cache: 'no-store'})
                .then(function (r) {
                    if (!r.ok) { throw new Error('bad status'); }
                    return {failures: 0, first: null, last: null};
                })
                .catch(function () {
                    return {failures: prev + 1, first: first || now, last: now};
                });
        }
        """,
        Output(STATE_ID, "data"),
        Input(PROBE_ID, "n_intervals"),
        State(STATE_ID, "data"),
    )

    app.clientside_callback(
        """
        function(state) {
            const n = (state && state.failures) || 0;
            const hidden = %s;
            const visible = %s;
            // A single missed beat is normal while a long synchronous callback holds the
            // worker. Only speak up once the server has really stopped answering.
            if (n < 2) { return [hidden, '', '']; }

            // ONE line, with a count — not one card per failed callback.
            const msg = 'Servidor nao respondeu \\u00b7 \\u00d7' + n + ' \\u00b7 ' + (state.last || '');
            const detail = 'Primeira ocorrencia: ' + (state.first || '?')
                + '\\nUltima: ' + (state.last || '?')
                + '\\nTentativas falhadas: ' + n
                + '\\n\\nMensagem (identica em todas): "the server did not respond".'
                + '\\nCausa provavel nesta maquina: o processo foi terminado por falta de'
                + ' memoria (OOM). Ve o painel de recursos no Setup.';
            return [visible, msg, detail];
        }
        """
        % (
            __import__("json").dumps(_HIDDEN),
            __import__("json").dumps(_VISIBLE),
        ),
        Output(BANNER_ID, "style"),
        Output(TEXT_ID, "children"),
        Output(DETAIL_ID, "children"),
        Input(STATE_ID, "data"),
    )
