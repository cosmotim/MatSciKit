"""
MatSciKit TUI — Interactive Terminal User Interface.

Launch with: matscikit (after pip install) or python -m matscikit_cosmotim.tui.app
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Center, Container, Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import (
    Button,
    DataTable,
    DirectoryTree,
    Footer,
    Header,
    Input,
    Label,
    ProgressBar,
    Select,
    Static,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PACKAGE_VERSION = "0.2.0"

SPLASH_ART = r"""
 __  __       _   ____       _ _  ___ _
|  \/  | __ _| |_/ ___|  ___(_) |/ (_) |_
| |\/| |/ _` | __\___ \ / __| | ' /| | __|
| |  | | (_| | |_ ___) | (__| | . \| | |_
|_|  |_|\__,_|\__|____/ \___|_|_|\_\_|\__|
"""

AUTHOR_INFO = """\
Author:  Yitian Wang
GitHub:  cosmotim
Email:   ywang1057@ucr.edu
"""

DATA_FORMATS = [
    ("PPMS HC  — Heat Capacity (.dat)", "ppms_hc"),
    ("PPMS TTO — Thermal Transport (.dat)", "ppms_tto"),
    ("DSC      — Differential Scanning Calorimetry (.csv)", "dsc"),
    ("LFA      — Laser Flash Analysis (.csv)", "lfa"),
]

# Column labels returned by each reader
FORMAT_COLUMNS = {
    "ppms_hc": ["Temperature (K)", "Cp (J/(g·K))", "Cp error"],
    "ppms_tto": ["Temperature (K)", "κ (W/m/K)", "κ error"],
    "dsc": ["Temperature (K)", "Cp (J/(g·K))"],
    "lfa": ["Temperature (K)", "Value", "Error"],
}


# ===================================================================
# Splash Screen
# ===================================================================


class SplashScreen(Screen):
    """Welcome screen shown on launch."""

    AUTO_DISMISS_SECONDS = 3.0

    CSS = """
    SplashScreen {
        align: center middle;
        background: $surface;
    }
    #splash-box {
        width: 64;
        height: auto;
        border: heavy $accent;
        padding: 1 2;
        background: $surface;
    }
    #splash-art {
        color: $accent;
        text-align: center;
    }
    #splash-version {
        color: $text-muted;
        text-align: center;
        margin-bottom: 1;
    }
    #splash-author {
        color: $text;
        text-align: center;
    }
    #splash-hint {
        color: $text-disabled;
        text-align: center;
        margin-top: 1;
    }
    """

    def compose(self) -> ComposeResult:
        with Center(), Container(id="splash-box"):
            yield Static(SPLASH_ART, id="splash-art")
            yield Static(f"v{PACKAGE_VERSION}", id="splash-version")
            yield Static(AUTHOR_INFO, id="splash-author")
            yield Static("Press any key to continue …", id="splash-hint")

    def on_mount(self) -> None:
        self.set_timer(self.AUTO_DISMISS_SECONDS, self._dismiss)

    def on_key(self) -> None:
        self._dismiss()

    def on_click(self) -> None:
        self._dismiss()

    def _dismiss(self) -> None:
        if self.is_current:
            self.app.pop_screen()


# ===================================================================
# Load Data Screen
# ===================================================================


class LoadDataScreen(Screen):
    """File browser + format selector + data preview."""

    BINDINGS = [
        Binding("escape", "go_back", "Back"),
    ]

    CSS = """
    LoadDataScreen {
        layout: horizontal;
        background: $surface;
    }
    #file-panel {
        width: 1fr;
        height: 100%;
        border-right: solid $accent;
        padding: 0 1;
    }
    #detail-panel {
        width: 1fr;
        height: 100%;
        padding: 0 1;
    }
    #format-select {
        margin: 1 0;
    }
    #load-btn {
        margin: 1 0;
        width: 100%;
    }
    #progress-bar {
        margin: 1 0;
    }
    #status-label {
        margin: 1 0;
        color: $accent;
    }
    #data-table {
        height: 1fr;
    }
    .section-title {
        color: $accent;
        text-style: bold;
        margin: 1 0 0 0;
    }
    DirectoryTree {
        height: 1fr;
    }
    """

    def __init__(self) -> None:
        super().__init__()
        self._selected_path: Path | None = None
        self._selected_format: str = "ppms_hc"

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        with Horizontal():
            with Vertical(id="file-panel"):
                yield Label("📂 Select Data File", classes="section-title")
                yield DirectoryTree(str(Path.home()), id="dir-tree")
            with Vertical(id="detail-panel"):
                yield Label("⚙️  Data Format", classes="section-title")
                yield Select(
                    [(label, value) for label, value in DATA_FORMATS],
                    value="ppms_hc",
                    id="format-select",
                )
                yield Button("Load Selected File", id="load-btn", variant="primary", disabled=True)
                yield ProgressBar(total=100, show_eta=False, id="progress-bar")
                yield Label("", id="status-label")
                yield Label("📊 Data Preview", classes="section-title")
                yield DataTable(id="data-table")
        yield Footer()

    @on(DirectoryTree.FileSelected, "#dir-tree")
    def file_selected(self, event: DirectoryTree.FileSelected) -> None:
        self._selected_path = Path(str(event.path))
        self.query_one("#load-btn", Button).disabled = False
        self.query_one("#status-label", Label).update(f"Selected: {self._selected_path.name}")

    @on(Select.Changed, "#format-select")
    def format_changed(self, event: Select.Changed) -> None:
        self._selected_format = str(event.value)

    @on(Button.Pressed, "#load-btn")
    def load_pressed(self) -> None:
        if self._selected_path:
            self._load_data()

    @work(thread=True)
    def _load_data(self) -> None:
        """Load data file using the appropriate reader."""
        progress = self.query_one("#progress-bar", ProgressBar)
        status = self.query_one("#status-label", Label)
        table = self.query_one("#data-table", DataTable)
        btn = self.query_one("#load-btn", Button)

        self.app.call_from_thread(btn.__setattr__, "disabled", True)
        self.app.call_from_thread(progress.update, progress=0)
        self.app.call_from_thread(status.update, f"Loading {self._selected_path.name} …")

        fmt = self._selected_format
        path = str(self._selected_path)

        # Simulate initial progress
        self.app.call_from_thread(progress.update, progress=20)

        try:
            # Import and call the correct reader
            if fmt == "ppms_hc":
                from MatSciKit_COSMOTIM.io import ppms_hc

                data = ppms_hc.read(path)
            elif fmt == "ppms_tto":
                from MatSciKit_COSMOTIM.io import ppms_tto

                data = ppms_tto.read(path)
            elif fmt == "dsc":
                from MatSciKit_COSMOTIM.io import dsc

                data = dsc.read(path)
            elif fmt == "lfa":
                from MatSciKit_COSMOTIM.io import lfa

                data = lfa.read(path)
            else:
                raise ValueError(f"Unknown format: {fmt}")

            self.app.call_from_thread(progress.update, progress=70)

            # Store loaded data on the app
            self.app.loaded_data = data  # type: ignore[attr-defined]
            self.app.loaded_format = fmt  # type: ignore[attr-defined]
            self.app.loaded_path = self._selected_path  # type: ignore[attr-defined]

            # Build summary
            n_rows, n_cols = data.shape
            t_min, t_max = float(np.nanmin(data[:, 0])), float(np.nanmax(data[:, 0]))
            summary = (
                f"✅ Loaded: {self._selected_path.name}\n"
                f"   Rows: {n_rows}  |  Columns: {n_cols}\n"
                f"   Temperature range: {t_min:.2f} – {t_max:.2f} K"
            )
            self.app.call_from_thread(status.update, summary)
            self.app.call_from_thread(progress.update, progress=90)

            # Populate data table
            columns = FORMAT_COLUMNS.get(fmt, [f"Col {i}" for i in range(n_cols)])

            def _populate_table():
                table.clear(columns=True)
                for col_name in columns[:n_cols]:
                    table.add_column(col_name)
                # Show first 50 rows
                for row in data[:50]:
                    table.add_row(*[f"{v:.6g}" for v in row])

            self.app.call_from_thread(_populate_table)
            self.app.call_from_thread(progress.update, progress=100)

        except Exception as exc:
            self.app.call_from_thread(status.update, f"❌ Error: {exc}")
            self.app.call_from_thread(progress.update, progress=0)

        self.app.call_from_thread(btn.__setattr__, "disabled", False)

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# About Screen
# ===================================================================


class AboutScreen(Screen):
    """Package information."""

    BINDINGS = [
        Binding("escape", "go_back", "Back"),
    ]

    CSS = """
    AboutScreen {
        align: center middle;
        background: $surface;
    }
    #about-box {
        width: 60;
        height: auto;
        border: heavy $accent;
        padding: 2 3;
    }
    #about-title {
        text-style: bold;
        color: $accent;
        text-align: center;
        margin-bottom: 1;
    }
    #about-text {
        text-align: center;
    }
    #about-back {
        margin-top: 2;
        width: 100%;
    }
    """

    def compose(self) -> ComposeResult:
        with Center(), Container(id="about-box"):
            yield Static("MatSciKit_COSMOTIM", id="about-title")
            yield Static(
                f"Version {PACKAGE_VERSION}\n\n"
                "A Python toolkit for materials science\n"
                "thermal property analysis.\n\n"
                "Pipelines:\n"
                "  1. Structure — crystal structure, Debye θ, sound velocities\n"
                "  2. Heat Capacity — Dulong-Petit, Debye model, Cp fitting\n"
                "  3. Thermal Conductivity — κ_min, MFP, LFA/TTO processing\n\n"
                f"{AUTHOR_INFO}",
                id="about-text",
            )
            yield Button("Back", id="about-back", variant="default")

    @on(Button.Pressed, "#about-back")
    def go_back_btn(self) -> None:
        self.app.pop_screen()

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# Main Menu Screen
# ===================================================================


class MainMenuScreen(Screen):
    """Primary navigation menu."""

    CSS = """
    MainMenuScreen {
        align: center middle;
        background: $surface;
    }
    #menu-box {
        width: 50;
        height: auto;
        border: heavy $accent;
        padding: 2 3;
    }
    #menu-title {
        text-style: bold;
        color: $accent;
        text-align: center;
        margin-bottom: 1;
    }
    .menu-btn {
        width: 100%;
        margin: 1 0;
    }
    """

    def compose(self) -> ComposeResult:
        with Center(), Container(id="menu-box"):
            yield Static("MatSciKit_COSMOTIM", id="menu-title")
            yield Button("📂  Load Data", id="btn-load", variant="primary", classes="menu-btn")
            yield Button(
                "🔬  Analysis Pipelines",
                id="btn-analysis",
                variant="default",
                classes="menu-btn",
                disabled=True,
            )
            yield Button("ℹ️   About", id="btn-about", variant="default", classes="menu-btn")
            yield Button("🚪  Quit", id="btn-quit", variant="error", classes="menu-btn")

    def on_screen_resume(self) -> None:
        """Enable analysis button if data has been loaded."""
        has_data = getattr(self.app, "loaded_data", None) is not None
        self.query_one("#btn-analysis", Button).disabled = not has_data

    @on(Button.Pressed, "#btn-load")
    def open_load(self) -> None:
        self.app.push_screen(LoadDataScreen())

    @on(Button.Pressed, "#btn-analysis")
    def open_analysis(self) -> None:
        self.app.push_screen(AnalysisScreen())

    @on(Button.Pressed, "#btn-about")
    def open_about(self) -> None:
        self.app.push_screen(AboutScreen())

    @on(Button.Pressed, "#btn-quit")
    def quit_app(self) -> None:
        self.app.exit()


# ===================================================================
# Analysis Screen (placeholder — will be expanded)
# ===================================================================


class AnalysisScreen(Screen):
    """Pipeline selection screen."""

    BINDINGS = [
        Binding("escape", "go_back", "Back"),
    ]

    CSS = """
    AnalysisScreen {
        align: center middle;
        background: $surface;
    }
    #analysis-box {
        width: 55;
        height: auto;
        border: heavy $accent;
        padding: 2 3;
    }
    #analysis-title {
        text-style: bold;
        color: $accent;
        text-align: center;
        margin-bottom: 1;
    }
    #analysis-info {
        color: $text-muted;
        text-align: center;
        margin-bottom: 1;
    }
    .analysis-btn {
        width: 100%;
        margin: 1 0;
    }
    """

    def compose(self) -> ComposeResult:
        with Center(), Container(id="analysis-box"):
            yield Static("🔬 Analysis Pipelines", id="analysis-title")
            yield Static("", id="analysis-info")
            yield Button(
                "1. Structure — Debye θ, Sound Velocities",
                id="btn-pipe1",
                variant="primary",
                classes="analysis-btn",
            )
            yield Button(
                "2. Heat Capacity — Dulong-Petit, Debye, Cp",
                id="btn-pipe2",
                variant="primary",
                classes="analysis-btn",
            )
            yield Button(
                "3. Thermal Conductivity — κ_min, MFP, LFA/TTO",
                id="btn-pipe3",
                variant="primary",
                classes="analysis-btn",
            )
            yield Button("← Back", id="btn-back", variant="default", classes="analysis-btn")

    def on_mount(self) -> None:
        loaded_path = getattr(self.app, "loaded_path", None)
        if loaded_path:
            self.query_one("#analysis-info", Static).update(f"Data loaded: {loaded_path.name}")

    @on(Button.Pressed, "#btn-back")
    def go_back_btn(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#btn-pipe1")
    def open_pipe1(self) -> None:
        self.app.push_screen(Pipeline1Screen())

    @on(Button.Pressed, "#btn-pipe2")
    def open_pipe2(self) -> None:
        self.app.push_screen(Pipeline2Screen())

    @on(Button.Pressed, "#btn-pipe3")
    def open_pipe3(self) -> None:
        self.app.push_screen(Pipeline3Screen())

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# Pipeline 1: Structure — Debye θ & Sound Velocities
# ===================================================================


class Pipeline1Screen(Screen):
    """Low-temperature Cp fitting for Debye temperature."""

    BINDINGS = [Binding("escape", "go_back", "Back")]

    CSS = """
    Pipeline1Screen {
        background: $surface;
    }
    .param-row {
        layout: horizontal;
        height: 3;
        margin: 0 1;
    }
    .param-label {
        width: 24;
        content-align: right middle;
        padding-right: 1;
        color: $text;
    }
    .param-input {
        width: 1fr;
    }
    #results-box {
        margin: 1 2;
        padding: 1 2;
        border: round $accent;
        height: auto;
    }
    #run-btn {
        margin: 1 2;
        width: auto;
    }
    """

    def compose(self) -> ComposeResult:

        yield Header(show_clock=True)
        yield Label("🔬 Pipeline 1: Debye Temperature from Low-T Cp", classes="section-title")
        with Horizontal(classes="param-row"):
            yield Label("n_density (atoms/m³):", classes="param-label")
            yield Input(placeholder="e.g. 7.63e28", id="inp-ndensity", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Density (kg/m³):", classes="param-label")
            yield Input(placeholder="e.g. 6870", id="inp-density", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("T min (K):", classes="param-label")
            yield Input(placeholder="3.0", id="inp-tmin", value="3.0", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("T max (K):", classes="param-label")
            yield Input(placeholder="10.0", id="inp-tmax", value="10.0", classes="param-input")
        yield Button("▶ Run Analysis", id="run-btn", variant="primary")
        yield Static("", id="results-box")
        yield Footer()

    @on(Button.Pressed, "#run-btn")
    def run_analysis(self) -> None:

        data = getattr(self.app, "loaded_data", None)
        getattr(self.app, "loaded_format", None)
        results = self.query_one("#results-box", Static)

        if data is None:
            results.update("❌ No data loaded. Load Cp data first.")
            return

        try:
            n_density = float(self.query_one("#inp-ndensity", Input).value)
            density = float(self.query_one("#inp-density", Input).value)
            t_min = float(self.query_one("#inp-tmin", Input).value)
            t_max = float(self.query_one("#inp-tmax", Input).value)
        except ValueError:
            results.update("❌ Please fill in all parameters with valid numbers.")
            return

        try:
            from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting

            T, Cp, Cp_err = data[:, 0], data[:, 1], data[:, 2]
            theta_D, v_s, theta_D_err, v_s_err = low_t_fitting.fit(
                T, Cp, Cp_err, n_density, density, t_range=(t_min, t_max)
            )
            results.update(
                "✅ Pipeline 1 Results\n"
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                f"  Debye temperature:  θ_D = {theta_D:.2f} ± {theta_D_err:.2f} K\n"
                f"  Sound velocity:     v_s = {v_s:.2f} ± {v_s_err:.2f} m/s\n"
                f"  Fitting range:      {t_min:.1f} – {t_max:.1f} K"
            )
        except Exception as exc:
            results.update(f"❌ Error: {exc}")

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# Pipeline 2: Heat Capacity — Dulong-Petit, Debye Model
# ===================================================================


class Pipeline2Screen(Screen):
    """Heat capacity analysis: Dulong-Petit limit and Debye model."""

    BINDINGS = [Binding("escape", "go_back", "Back")]

    CSS = """
    Pipeline2Screen {
        background: $surface;
    }
    .param-row {
        layout: horizontal;
        height: 3;
        margin: 0 1;
    }
    .param-label {
        width: 24;
        content-align: right middle;
        padding-right: 1;
        color: $text;
    }
    .param-input {
        width: 1fr;
    }
    #results-box {
        margin: 1 2;
        padding: 1 2;
        border: round $accent;
        height: auto;
    }
    #run-btn {
        margin: 1 2;
        width: auto;
    }
    """

    def compose(self) -> ComposeResult:

        yield Header(show_clock=True)
        yield Label("🔬 Pipeline 2: Heat Capacity Analysis", classes="section-title")
        with Horizontal(classes="param-row"):
            yield Label("n_density (atoms/m³):", classes="param-label")
            yield Input(placeholder="e.g. 7.63e28", id="inp-ndensity", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Density (kg/m³):", classes="param-label")
            yield Input(placeholder="e.g. 6870", id="inp-density", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Sound velocity (m/s):", classes="param-label")
            yield Input(placeholder="e.g. 2841 (optional)", id="inp-vs", classes="param-input")
        yield Button("▶ Run Analysis", id="run-btn", variant="primary")
        yield Static("", id="results-box")
        yield Footer()

    @on(Button.Pressed, "#run-btn")
    def run_analysis(self) -> None:

        results = self.query_one("#results-box", Static)

        try:
            n_density = float(self.query_one("#inp-ndensity", Input).value)
            density = float(self.query_one("#inp-density", Input).value)
        except ValueError:
            results.update("❌ Please provide n_density and density.")
            return

        try:
            from MatSciKit_COSMOTIM.heat_capacity import debye, dulong_petit

            dp = dulong_petit.calculate(n_density, density)

            lines = [
                "✅ Pipeline 2 Results",
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                f"  Dulong-Petit limit:  {dp:.4f} J/(g·K)",
            ]

            vs_str = self.query_one("#inp-vs", Input).value.strip()
            if vs_str:
                v_s = float(vs_str)
                theta_D = debye.from_velocity(v_s, n_density)
                lines.append(f"  Debye θ from v_s:    {theta_D:.2f} K")

            results.update("\n".join(lines))
        except Exception as exc:
            results.update(f"❌ Error: {exc}")

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# Pipeline 3: Thermal Conductivity — κ_min, MFP
# ===================================================================


class Pipeline3Screen(Screen):
    """Thermal conductivity analysis: κ_min and mean free path."""

    BINDINGS = [Binding("escape", "go_back", "Back")]

    CSS = """
    Pipeline3Screen {
        background: $surface;
    }
    .param-row {
        layout: horizontal;
        height: 3;
        margin: 0 1;
    }
    .param-label {
        width: 24;
        content-align: right middle;
        padding-right: 1;
        color: $text;
    }
    .param-input {
        width: 1fr;
    }
    #results-box {
        margin: 1 2;
        padding: 1 2;
        border: round $accent;
        height: auto;
    }
    #run-btn {
        margin: 1 2;
        width: auto;
    }
    """

    def compose(self) -> ComposeResult:

        yield Header(show_clock=True)
        yield Label(
            "🔬 Pipeline 3: Thermal Conductivity Analysis",
            classes="section-title",
        )
        with Horizontal(classes="param-row"):
            yield Label("n_density (atoms/m³):", classes="param-label")
            yield Input(placeholder="e.g. 7.63e28", id="inp-ndensity", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Debye θ (K):", classes="param-label")
            yield Input(placeholder="e.g. 312.4", id="inp-thetad", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Sound velocity (m/s):", classes="param-label")
            yield Input(placeholder="e.g. 2841", id="inp-vs", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("T for κ_min (K):", classes="param-label")
            yield Input(placeholder="e.g. 300", id="inp-temp", value="300", classes="param-input")
        with Horizontal(classes="param-row"):
            yield Label("Porosity (0-1):", classes="param-label")
            yield Input(
                placeholder="e.g. 0.03", id="inp-porosity", value="0.0", classes="param-input"
            )
        yield Button("▶ Run Analysis", id="run-btn", variant="primary")
        yield Static("", id="results-box")
        yield Footer()

    @on(Button.Pressed, "#run-btn")
    def run_analysis(self) -> None:

        results = self.query_one("#results-box", Static)

        try:
            n_density = float(self.query_one("#inp-ndensity", Input).value)
            theta_D = float(self.query_one("#inp-thetad", Input).value)
            v_s = float(self.query_one("#inp-vs", Input).value)
            T = float(self.query_one("#inp-temp", Input).value)
            porosity = float(self.query_one("#inp-porosity", Input).value)
        except ValueError:
            results.update("❌ Please fill in all parameters.")
            return

        try:
            from MatSciKit_COSMOTIM.thermal_conductivity import (
                cahill,
                mean_free_path,
                porosity_correction,
            )

            kappa_min = cahill.minimum_tc(T, n_density, theta_D, v_s)

            lines = [
                "✅ Pipeline 3 Results",
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                f"  κ_min at {T:.0f} K:    {kappa_min:.4f} W/(m·K)",
            ]

            if porosity > 0:
                kappa_corrected = porosity_correction.correct(kappa_min, porosity)
                lines.append(f"  κ (porosity={porosity:.2f}): {kappa_corrected:.4f} W/(m·K)")

            # If TTO data is loaded, compute MFP
            data = getattr(self.app, "loaded_data", None)
            fmt = getattr(self.app, "loaded_format", None)
            if data is not None and fmt == "ppms_tto":
                mfp = mean_free_path.calculate(data[:, 0], data[:, 1], theta_D, v_s)
                mfp_nm = mfp * 1e9
                lines.extend(
                    [
                        "",
                        "  Mean free path (from loaded TTO data):",
                        f"    Min: {np.nanmin(mfp_nm):.1f} nm  Max: {np.nanmax(mfp_nm):.1f} nm",
                    ]
                )

            results.update("\n".join(lines))
        except Exception as exc:
            results.update(f"❌ Error: {exc}")

    def action_go_back(self) -> None:
        self.app.pop_screen()


# ===================================================================
# Main App
# ===================================================================


class MatSciKitApp(App):
    """MatSciKit interactive terminal application."""

    TITLE = "MatSciKit"
    SUB_TITLE = f"v{PACKAGE_VERSION}"
    CSS = """
    Screen {
        background: $surface;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit", show=True),
    ]

    # State: set by LoadDataScreen when data is loaded
    loaded_data: np.ndarray | None = None
    loaded_format: str | None = None
    loaded_path: Path | None = None

    def on_mount(self) -> None:
        self.push_screen(MainMenuScreen())
        self.push_screen(SplashScreen())


def main() -> None:
    """Entry point for the matscikit command."""
    app = MatSciKitApp()
    app.run()


if __name__ == "__main__":
    main()
