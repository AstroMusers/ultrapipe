import ast
import csv
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest
from astropy.coordinates import SkyCoord


ROOT = Path(__file__).resolve().parent


def test_synthetic_field_notebook_uses_cost_pipeline_and_keeps_paths_private():
    path = ROOT / 'examples' / 'synthetic_field_ranking.ipynb'
    print('Reading from examples/synthetic_field_ranking.ipynb...')
    notebook = json.loads(path.read_text(encoding='utf-8'))
    code = '\n'.join(
        ''.join(cell['source'])
        for cell in notebook['cells']
        if cell['cell_type'] == 'code'
    )
    document = json.dumps(notebook)

    assert 'from utils import cost_function' in code
    assert 'score = cost_function(' in code
    assert 'save_bokeh_html(' in code
    assert 'get_visuals_path().parent / "examples" / "visuals"' in code
    assert 'explicitly synthetic' in document
    assert '/' + 'Users/' not in document
    static_path = ROOT / 'examples' / 'visuals' / 'ultrapipe_synthetic_field_scores.png'
    interactive_path = ROOT / 'examples' / 'visuals' / 'ultrapipe_synthetic_field_scores.html'
    assert static_path.is_file() and static_path.stat().st_size > 0
    print(f'Reading from {interactive_path.relative_to(ROOT)}...')
    interactive_html = interactive_path.read_text(encoding='utf-8')
    assert 'Synthetic ULTRASAT field scores' in interactive_html
    assert '/' + 'Users/' not in interactive_html


def test_utils_imports_without_online_model_dependencies(monkeypatch):
    for module_name in ("dustmaps", "dustmaps.bayestar", "expecto"):
        monkeypatch.setitem(sys.modules, module_name, None)
    spec = importlib.util.spec_from_file_location(
        "ultrapipe_utils_without_online_models", ROOT / "utils.py"
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)

    score = module.cost_function(
        0.0,
        0.0,
        7.0,
        [0.0],
        [0.0],
        [20.0],
        [9600.0],
        ROOT / "systematics" / "SatMag.csv",
        ROOT / "systematics" / "LimMag.csv",
        reddening=[0.0],
        target_weight_list=[1.0],
        data_weight_list=[1.0],
    )

    assert score == 0.25
    with pytest.raises(ImportError, match="dustmaps"):
        module.dust_map(resolution=1)
    with pytest.raises(ImportError, match="expecto"):
        module.TESS_ULTRASAT_color_transform([9600.0], [4.0])


@pytest.fixture
def cost_helpers():
    path = ROOT / 'utils.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    functions = {'find_nearest', 'SatMag', 'LimMag', 'cost_function'}
    constants = {'radial_pos', 'source_T'}
    selected = [node for node in tree.body if (
        isinstance(node, ast.FunctionDef) and node.name in functions
    ) or (
        isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in constants for target in node.targets)
    )]
    from tdpy.verbosity import print as verbosity_print

    namespace = {'np': np, 'csv': csv, 'math': math, 'SkyCoord': SkyCoord, 'print': verbosity_print}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def score_target(helpers, center, target, magnitude=20.0, reddening=0.0, target_weight=1.0, data_weight=1.0):
    return helpers['cost_function'](
        *center, 7.0, [target[0]], [target[1]], [magnitude], [9600.0],
        ROOT / 'systematics/SatMag.csv', ROOT / 'systematics/LimMag.csv',
        reddening=[reddening], target_weight_list=[target_weight], data_weight_list=[data_weight],
    )


@pytest.mark.parametrize(('center', 'target'), [
    ((1.0, 0.0), (359.0, 0.0)), ((359.0, 0.0), (1.0, 0.0)),
    ((0.0, 80.0), (10.0, 80.0)), ((0.0, -80.0), (10.0, -80.0)),
    ((0.0, 90.0), (180.0, 89.0)), ((0.0, -90.0), (180.0, -89.0)),
])
def test_field_membership_uses_spherical_wrap_safe_separation(cost_helpers, center, target):
    separation = SkyCoord(*center, unit='deg').separation(SkyCoord(*target, unit='deg')).degree  # [deg]
    assert separation < 7.0
    assert score_target(cost_helpers, center, target) == 0.25


@pytest.mark.parametrize('target', [(10.0, 0.0), (0.0, 10.0), (180.0, 0.0)])
def test_field_membership_retains_outside_target_sentinel(cost_helpers, target):
    assert score_target(cost_helpers, (0.0, 0.0), target) == -1


@pytest.mark.parametrize(('magnitude', 'expected'), [(13.0, 0.5), (15.0, 0.5), (20.0, 0.25), (24.0, -1)])
def test_valid_magnitude_bins_keep_their_original_scores(cost_helpers, magnitude, expected):
    assert score_target(cost_helpers, (0.0, 0.0), (0.0, 0.0), magnitude=magnitude) == expected


def test_field_score_keeps_weights_and_reddening(cost_helpers):
    assert score_target(cost_helpers, (0.0, 0.0), (0.0, 0.0), target_weight=2.0, data_weight=2.0) == 2.0
    assert score_target(cost_helpers, (0.0, 0.0), (0.0, 0.0), magnitude=15.0, reddening=1.0) == 0.25


def test_readme_uses_distinct_real_saturation_and_limiting_tables(cost_helpers):
    path = ROOT / 'README.rst'
    print(f'Reading from {path}...')
    readme = path.read_text()
    paths = {}
    for name, filename in (('SAT_Mag_dir', 'SatMag.csv'), ('LIM_Mag_dir', 'LimMag.csv')):
        assignment = f"{name} = get_repository_path() / 'systematics' / '{filename}'"
        assert assignment in readme
        paths[name] = ROOT / 'systematics' / filename
    assert paths['SAT_Mag_dir'] == ROOT / 'systematics/SatMag.csv'
    assert paths['LIM_Mag_dir'] == ROOT / 'systematics/LimMag.csv'
    saturation = cost_helpers['SatMag'](0.0, 9600.0, paths['SAT_Mag_dir'])
    limiting = cost_helpers['LimMag'](0.0, 9600.0, paths['LIM_Mag_dir'])
    assert saturation < 20.0 < limiting
    assert score_target(cost_helpers, (0.0, 0.0), (0.0, 0.0)) == 0.25


def test_field_score_keeps_empty_catalog_sentinel(cost_helpers):
    actual = cost_helpers['cost_function'](
        0.0, 0.0, 7.0, [], [], [], [],
        ROOT / 'systematics/SatMag.csv', ROOT / 'systematics/LimMag.csv',
        reddening=[], target_weight_list=[], data_weight_list=[],
    )
    assert actual == -1


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_sky_map_display_honors_theme_without_changing_values(monkeypatch, typeplotback, inherit):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import to_rgba
    from tdpy.plotting import plot_background_colors, plot_style

    path = ROOT / 'utils.py'
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == 'Plot_Cartesian_Sky_Map')
    namespace = {'plt': plt}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
    values = np.arange(12.0).reshape(3, 4)
    original_values = values.copy()
    background, foreground = plot_background_colors(typeplotback)
    displayed = []

    def inspect_display():
        figure = plt.gcf()
        displayed.append(figure)
        assert figure.get_facecolor() == to_rgba(background)
        for axis in figure.axes:
            assert axis.get_facecolor() == to_rgba(background)
            assert axis.xaxis.label.get_color() == foreground
            assert axis.yaxis.label.get_color() == foreground
        np.testing.assert_array_equal(figure.axes[0].images[0].get_array(), values)
        assert list(figure.axes[0].images[0].get_extent()) == [0, 360, 90, -90]
        FigureCanvasAgg(figure).draw()
        np.testing.assert_array_equal(np.asarray(figure.canvas.buffer_rgba())[0, 0], np.array(to_rgba(background)) * 255)

    monkeypatch.setattr(plt, 'show', inspect_display)
    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original_settings = dict(plt.rcParams)
        figure = plt.figure()
        try:
            assert namespace['Plot_Cartesian_Sky_Map'](values, 'Score', **options) is None
            assert dict(plt.rcParams) == original_settings
            assert len(displayed) == 1
            np.testing.assert_array_equal(values, original_values)
        finally:
            plt.close(figure)


def test_sky_map_saves_to_explicit_output_path(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    path = ROOT / 'utils.py'
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == 'Plot_Cartesian_Sky_Map')
    namespace = {'plt': plt}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
    monkeypatch.setattr(plt, 'show', lambda: None)
    output_path = tmp_path / 'visuals' / 'sky_map.png'

    assert namespace['Plot_Cartesian_Sky_Map'](np.arange(12.0).reshape(3, 4), output_path=output_path) is None
    assert output_path.is_file()
    assert output_path.stat().st_size > 1000
    plt.close('all')
