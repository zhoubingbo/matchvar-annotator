# 发布到 PyPI

本文说明如何把 `matchvar-annotator` 发布到 [TestPyPI](https://test.pypi.org/) 和 [PyPI](https://pypi.org/)。

## 1. 发布前准备

- 在 PyPI 和 TestPyPI 注册账号，并开启两步验证（2FA）。
- 确认包名 `matchvar-annotator` 可用，或你已经是该项目的维护者：
  - 若项目已存在，新版本号必须大于 PyPI 上的最新版本，否则上传会返回 `File already exists`。
- 版本号只有一个来源：`matchvar_annotator/__init__.py` 中的 `__version__`。
  `setup.py` 会自动读取它，不需要再改 `setup.py`。
- 需要构建工具（Python 3.8+）：

```bash
python -m pip install --upgrade pip setuptools wheel build twine
```

也可以直接安装开发依赖：

```bash
pip install -r requirements-dev.txt
```

## 2. 本地构建与检查

在项目根目录执行：

```bash
# 清理旧产物
rm -rf build dist *.egg-info

# 构建 sdist 和 wheel
python -m build

# 检查元数据、README 渲染和文件完整性
python -m twine check dist/*
```

`dist/` 下应得到类似文件：

```
dist/matchvar_annotator-1.2.0-py3-none-any.whl
dist/matchvar_annotator-1.2.0.tar.gz
```

建议先检查 wheel 内容，确认没有把 `tests/`、缓存和本机数据打进去：

```bash
python -m zipfile -l dist/matchvar_annotator-1.2.0-py3-none-any.whl
```

## 3. 上传到 TestPyPI（推荐先做）

先在 [TestPyPI](https://test.pypi.org/manage/account/token/) 创建 API token，然后执行：

```bash
python -m twine upload --repository testpypi dist/*
```

用户名固定为 `__token__`，密码填 `pypi-` 开头的 token。也可以用环境变量：

```bash
export TWINE_USERNAME=__token__
export TWINE_PASSWORD=pypi-xxxxxxxxxxxxxxxx
python -m twine upload --repository testpypi dist/*
```

验证安装（TestPyPI 依赖不全，需加官方源）：

```bash
python -m venv /tmp/mv-test
source /tmp/mv-test/bin/activate
pip install --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple \
  matchvar-annotator==1.2.0
matchvar-annotator --help
```

## 4. 上传到正式 PyPI

在 [PyPI](https://pypi.org/manage/account/token/) 创建 API token 后执行：

```bash
python -m twine upload dist/*
```

同样使用 `__token__` 作为用户名。发布成功后检查：

- <https://pypi.org/project/matchvar-annotator/>
- 确认版本号、README、License、项目链接、可选依赖 `[bb]` 均显示正确。

### 使用 `~/.pypirc`（可选）

复制示例文件并填入 token：

```bash
cp .pypirc.example ~/.pypirc
chmod 600 ~/.pypirc
```

之后可直接使用：

```bash
python -m twine upload --repository pypi dist/*
python -m twine upload --repository testpypi dist/*
```

不要把 `.pypirc` 或 token 提交到 git（仓库的 `.gitignore` 已忽略 `.pypirc`）。

## 5. 通过 GitHub Actions 自动发布（推荐）

仓库已包含 `.github/workflows/publish.yml`，使用 PyPI 的 **Trusted Publishing（OIDC）**，无需在仓库里保存 token。

### 5.1 在 PyPI 配置发布者

在 PyPI 项目页 → **Publishing** → **Add a new publisher**（新项目用 **Pending publisher**）填写：

| 字段 | 值 |
|------|-----|
| PyPI project name | `matchvar-annotator` |
| Owner | `zhoubingbo` |
| Repository name | `matchvar-annotator` |
| Workflow name | `publish.yml` |
| Environment name | `pypi` |

在 TestPyPI 做同样配置，Environment name 填 `testpypi`，仓库 URL 使用 <https://test.pypi.org/legacy/>。

### 5.2 在 GitHub 创建 Environment

仓库 **Settings → Environments** 新建：

- `pypi`
- `testpypi`

可选：为 `pypi` 添加 Required reviewers，发布前需要人工批准。

### 5.3 触发发布

- **正式发布**：在 GitHub 创建 Release（建议 tag 形如 `v1.2.0`），工作流会自动构建并上传 PyPI。
- **测试发布**：Actions → **Publish to PyPI** → **Run workflow**，目标选择 `testpypi`。

## 6. 版本发布清单

1. 修改 `matchvar_annotator/__init__.py` 中的 `__version__`。
2. 在 `CHANGELOG.md` 记录变更。
3. 本地运行测试：`pytest tests/`。
4. 本地构建并检查：`python -m build && python -m twine check dist/*`。
5. 上传 TestPyPI 并验证安装。
6. 打 tag 并创建 GitHub Release，等待自动发布到 PyPI。
7. 验证 `pip install matchvar-annotator` 和 `matchvar-annotator --help`。

## 7. 常见问题

| 错误 | 原因与处理 |
|------|------------|
| `HTTPError: 403 Forbidden` | token 无该项目权限、用户名不是 `__token__`，或 2FA/token 作用域不正确 |
| `File already exists` | 同版本文件已上传；PyPI 不允许覆盖，需提升 `__version__` 后重新构建 |
| `InvalidDistribution: metadata is missing required fields` | `setup.py` 元数据不完整，运行 `python -m twine check dist/*` 查看详情 |
| `README` 在 PyPI 上排版错乱 | 确认 `long_description_content_type="text/markdown"`，并用 `twine check` 检查 |
| 构建时报找不到 `setuptools>=64` | 升级构建环境：`python -m pip install -U setuptools wheel build`；国内镜像可加 `PIP_INDEX_URL=https://pypi.org/simple` |
| wheel 里包含 `tests/` | 确认 `setup.py` 使用 `find_packages(exclude=("tests", "tests.*"))` |
| `matchvar-annotator` 名称被占用 | 需在 PyPI 申请项目所有权，或改用其他包名并同步更新所有文档与脚本 |
