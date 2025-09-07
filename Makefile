# MATCHVAR注释器 Makefile

.PHONY: help install install-dev test lint format clean build dist upload docs

help:  ## 显示帮助信息
	@echo "可用的命令:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## 安装包
	pip install -e .

install-dev:  ## 安装开发依赖
	pip install -e ".[dev]"

test:  ## 运行测试
	pytest tests/ -v --cov=matchvar_annotator --cov-report=html --cov-report=term

test-fast:  ## 快速测试
	pytest tests/ -v -x

lint:  ## 代码检查
	flake8 matchvar_annotator/ tests/
	mypy matchvar_annotator/

format:  ## 代码格式化
	black matchvar_annotator/ tests/
	isort matchvar_annotator/ tests/

format-check:  ## 检查代码格式
	black --check matchvar_annotator/ tests/
	isort --check-only matchvar_annotator/ tests/

clean:  ## 清理临时文件
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .mypy_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean  ## 构建包
	python -m build

dist: build  ## 创建分发包
	@echo "分发包已创建在 dist/ 目录中"

upload-test: dist  ## 上传到测试PyPI
	twine upload --repository testpypi dist/*

upload: dist  ## 上传到PyPI
	twine upload dist/*

docs:  ## 构建文档
	cd docs && make html

docs-serve: docs  ## 本地服务文档
	cd docs/_build/html && python -m http.server 8000

check: format-check lint test  ## 运行所有检查

ci:  ## CI/CD检查
	pip install -e ".[dev]"
	black --check matchvar_annotator/ tests/
	flake8 matchvar_annotator/ tests/
	mypy matchvar_annotator/
	pytest tests/ -v --cov=matchvar_annotator --cov-report=xml

version:  ## 显示版本信息
	@python -c "import matchvar_annotator; print(f'版本: {matchvar_annotator.__version__}')"

info:  ## 显示包信息
	@python -c "import matchvar_annotator; print(f'包名: {matchvar_annotator.__name__}'); print(f'版本: {matchvar_annotator.__version__}'); print(f'作者: {matchvar_annotator.__author__}'); print(f'描述: {matchvar_annotator.__description__}')"
